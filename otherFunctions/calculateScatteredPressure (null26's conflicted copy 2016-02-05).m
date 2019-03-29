function F_k = calculateScatteredPressure(varCol, U, k, P_far, noWavesVec)

n_res = 12;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;


% Check orientation of NURBS object (assuming the object is orientable)
rightHandedOrientation = findOrientation(varCol);

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

noElemsXi = length(unique(varCol.nurbs.knots{1})) - 1;
noElemsEta = length(unique(varCol.nurbs.knots{2})) - 1;

noWavesInXiDir = noWavesVec(1);
noWavesInEtaDir = noWavesVec(2);

noGpXi = ceil(n_res*noWavesInXiDir/noElemsXi);
noGpEta = ceil(n_res*noWavesInEtaDir/noElemsEta);

if noGpXi <= p+1
    noGpXi = p+1;
end
if noGpEta <= q+1
    noGpEta = q+1;
end

[W2D,Q2D] = gaussianQuadNURBS(noGpXi,noGpEta);  

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

        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);

        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        crossProd = cross(J(:,1),J(:,2));
        if rightHandedOrientation
            normal = crossProd/norm(crossProd);
        else
            normal = -crossProd/norm(crossProd);
        end
        
        v = evaluateNURBS(nurbs, [xi, eta, zeta]);           
        
        diffr = zeros(size(P_far));
        diffr(1,:,:) = v(1)-P_far(1,:,:);
        diffr(2,:,:) = v(2)-P_far(2,:,:);
        diffr(3,:,:) = v(3)-P_far(3,:,:);
        
        norm_rm_r = norm2(diffr);    
        greensFunc = exp(1i*k*norm_rm_r)./(4*pi*norm_rm_r);    
        dgreensFunc = greensFunc.*(1i*k*norm_rm_r - 1)./(norm_rm_r.^2) .*dot3(diffr, normal);

%         p_h = R_fun*U(sctr,:,:);
%         dp_h = dot3(matrix3Dprod(J'\[dRdxi; dRdeta; dRdzeta],U(sctr,:,:)),normal);
        p_h = R_fun*U(sctr,:);
        dp_h = dot3(J'\[dRdxi; dRdeta; dRdzeta]*U(sctr,:),normal);

        F_k = F_k + (p_h.*dgreensFunc - dp_h.*greensFunc)* norm(crossProd) * J_2 * wt; 
        
    end
%     if length(P_far) > 3
%         disp(['Completed ' num2str(i) ' out of ' num2str(length(innerSurfaceElements)) ' Elapsed time: ' num2str(toc)])
%     end
end