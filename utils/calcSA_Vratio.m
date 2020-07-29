function ratio = calcSA_Vratio(varCol)
error('Clean up this routine')

volume = measureNURBS(varCol.nurbs);


% Find elements on the inner surface for evaluation of
% backscattered pressure in far field
innerSurfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(1) == 0
        innerSurfaceElements = [innerSurfaceElements e];
    end
end

[W2D,Q2D] = gaussianQuadNURBS(p_xi+3,p_eta+3);  

area = 0;

parfor i = 1:length(innerSurfaceElements) %[8 7 3 4 5 6 2 1]% 
% for i = 1:length(innerSurfaceElements)
    e = innerSurfaceElements(i);
    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = 0;
        
        [~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, degree, Xi, Eta, Zeta, weights);

        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        crossProd = cross(J(:,1),J(:,2));
        J_1 = norm(crossProd);
        
        area = area + J_1 * J_2 * wt;
    end
end

ratio = area/volume;
