function Error = calcError1D(varCol, U)

%% Extract all needed data from options and varCol
Xi = varCol.patches{1}.nurbs.knots;
p_xi = varCol.patches{1}.nurbs.degree;
weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
element = varCol.element;

analytic = varCol.analytic;

%% Preallocation and initiallizations
[W1D,Q1D] = gaussianQuadNURBS(p_xi+10); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
Error = 0;
normalization = 0;
% parfor e = 1:noElems
for e = 1:noElems
    idXi = index(e);
    
    Xi_e = elRangeXi(idXi,:);
    
    J_2 = 0.5*(Xi_e(2)-Xi_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr);
    U_sctr = U(sctr);
    
    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi = parent2ParametricSpace(Xi_e,  pt(1));
        
        [R_fun, dRdxi] = NURBS1DBasis(xi, p_xi, Xi, weights);
        
        J = dRdxi*pts;
        J_1 = J;
%         dRdX = dRdxi/J_1;
        
        u_h = R_fun*U_sctr;
        X = R_fun*pts;
        u = analytic(X);
        Error = Error + abs(u-u_h)^2 * abs(J_1) * J_2 * wt;
        normalization = normalization + abs(u)^2 * abs(J_1) * J_2 * wt;
    end
end

Error = 100*sqrt(Error/normalization);
