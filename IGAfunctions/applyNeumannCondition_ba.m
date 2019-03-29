function F = applyNeumannCondition_ba(varCol,analytic, gradient,v_arr)

%% Extract all needed data from varCol
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
l = varCol.nurbs.number(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;

noCtrlPts = varCol.noCtrlPts;
noDofs = varCol.noDofs;

%% Preallocation and initiallizations
n_en = (p+1)*(q+1)*(r+1);

F_indices = zeros(n_en,noElems); 
Fvalues = zeros(n_en,noElems); 

[W3D,Q3D] = gaussianQuadNURBS(p+1+floor(16*p/n),q+1+floor(8*q/m),r+1+floor(4*r/l)); 
% [W3D,Q3D] = gaussianQuadNURBS(p+1,q+1,r+1); 

analytic_values = analytic(v_arr);
gradient_values = gradient(v_arr);
analytic_values = reshape(analytic_values, size(W3D,1), noElems, 1);
gradient_values = reshape(gradient_values, size(W3D,1), noElems, 3);
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    Zeta_e = elRangeZeta(idZeta,:);
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    
    f_e = zeros(n_en,1);
    
    for gp = 1:size(W3D,1)
        pt = Q3D(gp,:);
        wt = W3D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = parent2ParametricSpace(Zeta_e,pt(3));
        
        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
        
        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta; dRdzeta];
        
        p_analytic = analytic_values(gp,e);
        gradient_analytic = zeros(3,1);
        gradient_analytic(:) = gradient_values(gp,e,:);
        f_e = f_e + (dot2(dRdX, gradient_analytic).' + p_analytic*R_fun') * abs(J_1) * J_2 * wt;        
    end

    F_indices(:,e) = sctr';
    Fvalues(:,e) = f_e;
end
F = vectorAssembly(Fvalues,F_indices,noDofs);



