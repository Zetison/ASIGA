function [A, FF] = buildBEMmatrix_galerkin2(varCol)
error('Used buildCBEMmatrix instead')

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;

noElems = varCol.noElems;
index = varCol.index;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);



noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

P_inc = varCol.P_inc;

%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1);

Fvalues = zeros(n_en,noElems); 
Findices = zeros(n_en,noElems); 
spIdxRow = zeros(n_en^2,noElems);
spIdxCol = zeros(n_en^2,noElems);
Avalues = zeros(n_en^2,noElems); 
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    F_e = zeros(n_en,1);
    A_e = zeros(n_en);
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

        J = pts'*[dRdxi' dRdeta'];
        crossProd = cross(J(:,1),J(:,2));
        J_1 = norm(crossProd);

        x = R*pts;
        fact_x = J_1*J_2*wt;
        
        A_e = A_e - R.'*R*fact_x;
        F_e = F_e - R.'*P_inc(x).'*fact_x;
    end
    Avalues(:,e) = reshape(A_e, n_en^2, 1);
    
    Fvalues(:,e) = F_e;
    Findices(:,e) = sctr;
    spIdxRow(:,e) = copyVector(sctr,n_en,1);
    spIdxCol(:,e) = copyVector(sctr,n_en,2);
end
A = sparse(spIdxRow,spIdxCol,Avalues,noDofs,noDofs);
FF = vectorAssembly(Fvalues, Findices, noDofs);


