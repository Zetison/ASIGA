function [L2Error,H1error,H1sError] = calcError2D(varCol, U)

%% Extract all needed data from options and varCol
Xi = varCol.patches{1}.nurbs.knots{1};
Eta = varCol.patches{1}.nurbs.knots{2};
p_xi = varCol.patches{1}.nurbs.degree(1);
p_eta = varCol.patches{1}.nurbs.degree(2);
weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;

index = varCol.patches{1}.index;
noElems = varCol.patches{1}.noElems;
elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};
element = varCol.patches{1}.element;

analytic = varCol.analytic;
danalytic = varCol.danalytic;

%% Preallocation and initiallizations
[W2D,Q2D] = gaussianQuadNURBS(p_xi+5,p_eta+5); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
L2Error = 0;
H1error = 0;
H1sError = 0;
normalizationL2 = 0;
normalizationH1 = 0;
normalizationH1s = 0;
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    U_sctr = U(sctr);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi = parent2ParametricSpace(Xi_e,  pt(1));
        eta = parent2ParametricSpace(Eta_e,  pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
        
        J = pts'*[dRdxi' dRdeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta];
        
        u_h = R*U_sctr;
        du_h = dRdX*U_sctr;
        y = R*pts;
        u = analytic(y(1),y(2));
        du = danalytic(y(1),y(2));
        L2Error = L2Error + abs(u-u_h)^2 * abs(J_1) * J_2 * wt;
        H1error = H1error + (norm(du-du_h)^2 + abs(u-u_h)^2) * abs(J_1) * J_2 * wt;
        H1sError = H1sError + norm(du-du_h)^2 * abs(J_1) * J_2 * wt;
        normalizationL2 = normalizationL2 + abs(u)^2 * abs(J_1) * J_2 * wt;
        normalizationH1 = normalizationH1 + (norm(du)^2 + abs(u)^2) * abs(J_1) * J_2 * wt;
        normalizationH1s = normalizationH1s + norm(du)^2 * abs(J_1) * J_2 * wt;
    end
end

L2Error = 100*sqrt(L2Error/normalizationL2);
H1error = 100*sqrt(H1error/normalizationH1);
H1sError = 100*sqrt(H1sError/normalizationH1s);
