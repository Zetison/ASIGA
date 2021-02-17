function p_h = calculateScatteredPressure(varCol, P_far, useExtraQuadPts, computeFarField)

d_p = varCol{1}.patches{1}.nurbs.d_p;
d = varCol{1}.patches{1}.nurbs.d;
U = varCol{1}.U;
farFieldNormalPressFromSolid = varCol{1}.farFieldNormalPressFromSolid;
if numel(varCol) > 1 && (d_p == 2 || farFieldNormalPressFromSolid)
    [nodesSolid, ~, elementSolid] = meshBoundary(varCol{2},1);
    noDofs = varCol{2}.noDofs;
    Ux = varCol{2}.U(1:d:noDofs,:);
    Uy = varCol{2}.U(2:d:noDofs,:);
    Uz = varCol{2}.U(3:d:noDofs,:);
else
    nodesSolid = NaN;
    elementSolid = NaN;
    Ux = NaN;
    Uy = NaN;
    Uz = NaN;
end
knotVecs = varCol{1}.knotVecs;
weights = varCol{1}.weights;
controlPts = varCol{1}.controlPts;
elRange = varCol{1}.elRange;
noDofs = varCol{1}.noDofs;
rho = varCol{1}.rho;
if farFieldNormalPressFromSolid && d_p == 3
    degree = varCol{1}.degree(1:2); % assume p_xi is equal in all patches
    [zeta0Nodes, noElems, element, element2, index, pIndex] = meshBoundary(varCol{1},0);
else
    degree = varCol{1}.degree; % assume p_xi is equal in all patches
    index = varCol{1}.index;
    noElems = varCol{1}.noElems;
    element = varCol{1}.element;
    element2 = varCol{1}.element2;
    pIndex = varCol{1}.pIndex;
    zeta0Nodes = 1:noDofs;
end

BC = varCol{1}.BC;
k = varCol{1}.k;
omega = varCol{1}.omega;
method = varCol{1}.method;
if strcmp(method,'IENSG') && ~farFieldNormalPressFromSolid
    error('Not implemented')
end
dp_inc = varCol{1}.dp_inc_;

Phi_k = varCol{1}.Phi_k;
if d_p == 3 && ~farFieldNormalPressFromSolid
    surfaceElements = [];
    for e = 1:noElems
        idZeta = index(e,3);
        Zeta_e = elRange{3}(idZeta,:); % [zeta_k,zeta_k+1]                    
        if Zeta_e(1) == 0
            surfaceElements = [surfaceElements e];
        end
    end
else
    surfaceElements = 1:noElems;
end
solveForPtot = varCol{1}.solveForPtot;
exteriorSHBC = (strcmp(BC, 'SHBC') || strcmp(BC, 'NBC')) && numel(varCol) == 1;

extraGP = varCol{1}.extraGP;
if useExtraQuadPts
    noGp = degree+1+5;
else
    noGp = degree+3+extraGP(1:numel(degree));
end
[Q, W] = gaussTensorQuad(noGp); 
if numel(k) > 1 && size(P_far,1) > 1
    error('not implemented')
end
p_h = zeros(max([size(P_far,1),numel(k)]),1);

% for i = 1:length(surfaceElements) %
parfor i = 1:length(surfaceElements)
    e = surfaceElements(i);
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(2,2);
    for ii = 1:2
        Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
    end
    
    sctr = zeta0Nodes(element(e,:));
    pts = controlPts(sctr,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:);

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;

    xi = [parent2ParametricSpace(Xi_e, Q), zeros(size(Q,1),1)];
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1,crossProd] = getJacobian(R,pts,2);
    normals = crossProd./repmat(J_1,1,3);

    Y = R{1}*pts;  
    U_sctr = U(sctr,:);
    p_h_gp = R{1}*U_sctr;

    if exteriorSHBC
        dp_h_gp = -dp_inc(Y,normals); 
    else
        if d_p == 2 || farFieldNormalPressFromSolid
            sctrSolid = nodesSolid(elementSolid(e,:));
            dp_h_gp = rho*omega.^2.*( normals(:,1).*(R{1}*Ux(sctrSolid,:)) ...
                                     +normals(:,2).*(R{1}*Uy(sctrSolid,:)) ...
                                     +normals(:,3).*(R{1}*Uz(sctrSolid,:)));
            if ~solveForPtot
                dp_h_gp = dp_h_gp - dp_inc(Y,normals);
            end
        else
            dp_h_gp = zeros(size(p_h_gp));
            for gp = 1:size(W,1)
                dXdxi = R{2}*pts;
                dXdeta = R{3}*pts;
                dXdzeta = R{4}*pts;
                J = [dXdxi(gp,:).' dXdeta(gp,:).' dXdzeta(gp,:).'];
                dp_h_gp(gp,:) = (J'\[R{2}(gp,:); R{3}(gp,:); R{4}(gp,:)]*U_sctr).'*normals(gp,:).';
            end
        end
    end

    X = P_far./repmat(norm2(P_far),1,size(P_far,2));
    
    if computeFarField
        x_d_n = normals*X.';
        x_d_y = Y*X.';
        if solveForPtot
            p_h = p_h + (1i*k.*p_h_gp.*x_d_n.*exp(-1i*k.*x_d_y)).'* (J_1 * J_2 .* W);  
        else
            p_h = p_h + ((1i*k.*p_h_gp.*x_d_n + dp_h_gp).*exp(-1i*k.*x_d_y)).'* (J_1 * J_2 .* W);  
        end
    else
        xmy = reshape(P_far,[1,size(P_far,1),3]) - reshape(Y,[size(Y,1),1,3]);
        r = norm2(xmy);
        dPhidny =  Phi_k(r)./r.^2.*(1 - 1i*k.*r).*sum(xmy.*repmat(reshape(normals,size(normals,1),1,3),1,size(P_far,1),1),3);
        if solveForPtot
            p_h = p_h + (p_h_gp.*dPhidny).'* (J_1 * J_2 .* W); 
        else 
            p_h = p_h + (p_h_gp.*dPhidny - dp_h_gp.*Phi_k(r)).'* (J_1 * J_2 .* W);
        end 
    end
end
if numel(k) > 1
    p_h = p_h.';
end
if computeFarField
    p_h = -1/(4*pi)*p_h;
end



