function F = applyNeumannPulsatingSphere(varCol,k,P_0,R_o)

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

noElems = varCol.noElems;
index = varCol.index;
element = varCol.element;


noDofs_tot = varCol.noDofs_tot;

weights = varCol.weights;
controlPts = varCol.controlPts;


% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs_tot,1);        % external force vector

n_en = (p+1)*(q+1);
Fvalues = zeros(n_en,noElems);
indices = zeros(n_en,noElems);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
% parfor e = 1:noElems
for e = 1:noElems
    idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);          %  element scatter vector

    n_en = length(sctr);
    pts = controlPts(sctr,:);
    f_e = zeros(n_en,1);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2)); 
       
        deriv = -P_0*(1i*k-1/R_o);
        
        f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctr';
    Fvalues(:,e) = f_e;
end
F = vectorAssembly(Fvalues,indices,noDofs_tot);
