function F_k = calculateBackscatteredPressure_2D(varCol, U, noGpXi, k, P_far)

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;


% Check orientation of NURBS object (assuming the object is orientable)
rightHandedOrientation = findOrientation_2D(varCol);

% Find elements on the inner surface for evaluation of
% backscattered pressure in far field
innerSurfaceElements = [];
for e = 1:noElems
    idEta = index(e,2);
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]                    
    if Eta_e(1) == 0
        innerSurfaceElements = [innerSurfaceElements e];
    end
end

[W1D,Q1D] = gaussianQuadNURBS(noGpXi);  

if size(P_far,2) < size(U,2)
    error('The number of angles in P_far must be greater or equal to the number of angles in U')
end
if size(U,3) > 1
    error('This is not implemented, but ready to use if changed to matrix3Dprod(A,B)...')
end
F_k = zeros(1, size(P_far,2), size(P_far,3));

% for i = 1:length(innerSurfaceElements) %
parfor i = 1:length(innerSurfaceElements) %[8 7 3 4 5 6 2 1]% 
    tic
    e = innerSurfaceElements(i);
    idXi = index(e,1);

    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]

    J_2 = 0.5*(Xi_e(2)-Xi_e(1));

    sctr = element(e,:);

    pts = controlPts(sctr,:);

    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = 0;

        [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);

        J = pts'*[dRdxi' dRdeta'];
        crossProd = J(:,1);
        
        tangent = crossProd/norm(crossProd);
        normal = -[-tangent(2), tangent(1)]';
        
        v = evaluateNURBS(nurbs, [xi, eta]);           
        
        diffr = zeros(size(P_far));
        diffr(1,:,:) = v(1)-P_far(1,:);
        diffr(2,:,:) = v(2)-P_far(2,:);
        
        norm_rm_r = norm2(diffr);    
        greensFunc = 1i/4*cylindricalHankel1(0,k*norm_rm_r);    
%         dgreensFunc = 1i*k/4*dot3(diffr, normal).*cylindricalHankel1Deriv(0,k*norm_rm_r)./norm_rm_r;
        dgreensFunc = -1i*k/4*dot3(diffr, normal).*cylindricalHankel1(1,k*norm_rm_r)./norm_rm_r;
        
%         p_h = R_fun*U(sctr,:,:);
%         dp_h = dot3(matrix3Dprod(J'\[dRdxi; dRdeta; dRdzeta],U(sctr,:,:)),normal);
        p_h = R_fun*U(sctr,:);
        dp_h = dot3(J'\[dRdxi; dRdeta]*U(sctr,:),normal);
        

        F_k = F_k + (p_h.*dgreensFunc - dp_h.*greensFunc)* norm(crossProd) * J_2 * wt; 
    end
end