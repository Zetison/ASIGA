function p_h = calculateScatteredPressureKDT(varCol, P_far, plotFarField)


if ~strcmp(varCol.BC, 'SHBC')
    error('This is not implemented')
end
p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;

d_vec = varCol.d_vec;
k = varCol.k;

p_inc = varCol.p_inc;

extraGP = varCol.extraGP;

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);  
% [W2D,Q2D] = gaussianQuadNURBS(20,20);

p_h = zeros(size(P_far,1),length(k));
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:)); % New   

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));

        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

        J = pts'*[dRdxi' dRdeta'];
        crossProd = cross(J(:,1),J(:,2));
        J_1 = norm(crossProd);
        n = crossProd/J_1;

        Y = R*pts;

        p_h_gp = 2*p_inc(Y);
        p_h_gp(n.'*d_vec > 0) = 0;
        
        if plotFarField
            x_d_n = P_far*n./norm2(P_far);
            x_d_y = P_far*Y.'./norm2(P_far);
            p_h = p_h - 1/(4*pi)*1i*k.*p_h_gp.'.*x_d_n.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
        else
            xmy = P_far - repmat(Y,size(P_far,1),1);
            r = norm2(xmy);
            p_h = p_h + p_h_gp.'.*dPhi_kdny(xmy,r,n,k)* J_1 * J_2 * wt;  
        end
    end
end


