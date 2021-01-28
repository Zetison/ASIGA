function p_h = calculateScatteredPressure(varCol, Uc, P_far, useExtraQuadPts, computeFarField)

d_p = varCol{1}.patches{1}.nurbs.d_p;
d = varCol{1}.patches{1}.nurbs.d;
U = Uc{1};
if numel(varCol) > 1 && d_p == 2
    [~, ~, elementSolid] = meshBoundary(varCol{2},1);
    noDofs = varCol{2}.noDofs;
    Ux = Uc{2}(1:d:noDofs,:);
    Uy = Uc{2}(2:d:noDofs,:);
    Uz = Uc{2}(3:d:noDofs,:);
else
    elementSolid = NaN;
    Ux = NaN;
    Uy = NaN;
    Uz = NaN;
end
rho = varCol{1}.rho;
degree = varCol{1}.degree; % assume p_xi is equal in all patches
noDofs = varCol{1}.noDofs;
index = varCol{1}.index;
noElems = varCol{1}.noElems;
elRange = varCol{1}.elRange;
element = varCol{1}.element;
element2 = varCol{1}.element2;
weights = varCol{1}.weights;
controlPts = varCol{1}.controlPts;
knotVecs = varCol{1}.knotVecs;
pIndex = varCol{1}.pIndex;
BC = varCol{1}.BC;
k = varCol{1}.k;
omega = varCol{1}.omega;
method = varCol{1}.method;
if strcmp(method,'IENSG')
    A_2 = varCol{1}.A_2;
    x_0 = varCol{1}.x_0;
    Upsilon = varCol{1}.Upsilon;
    D = varCol{1}.D;
    N = varCol{1}.N;
else
    A_2 = NaN;
    x_0 = NaN;
    Upsilon = NaN;
    D = NaN;
    N = NaN;
end

dp_inc = varCol{1}.dp_inc;

Phi_k = varCol{1}.Phi_k;
if d_p == 3
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
if exteriorSHBC
    if solveForPtot
        homNeumanCond = true;
        dpdn = @(x,n) 0;
    else
        homNeumanCond = false;
        dpdn = @(x,n) -varCol{1}.dp_inc(x,n);
    end
else
    homNeumanCond = false;
    dpdn = varCol{1}.dpdn;
end

extraGP = varCol{1}.extraGP;
if useExtraQuadPts
    noGpXi = degree(1)+1+5;
    noGpEta = degree(2)+1+5;
else
    noGpXi = degree(1)+3+extraGP;
    noGpEta = degree(2)+3+extraGP;
end
[Q, W] = gaussTensorQuad([noGpXi,noGpEta]); 

p_h = zeros(size(P_far,1),numel(k));

% for i = 1:length(surfaceElements) %
parfor i = 1:length(surfaceElements)
    if d_p == 3
        e = surfaceElements(i);
    else
        e = i;
    end
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(2,2);
    for ii = 1:2
        Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;

    xi = [parent2ParametricSpace(Xi_e, Q), zeros(size(Q,1),1)];
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1,crossProd] = getJacobian(R,pts,2);
    normals = crossProd./repmat(J_1,1,3);

    Y = R{1}*pts;  
    if strcmp(method,'IENSG')
        Yt = (Y-x_0)*A_2.';

        r_a = evaluateProlateCoords(Yt,Upsilon);
        p_h_gp = zeros(size(W,1),size(U,2));
        for m = 1:N
            temp = 0;
            temp2 = 0;
            for mt = 1:N
                temp = temp + (1i*k - mt/r_a)*D(m,mt);
                temp2 = temp2 + D(m,mt);
            end
            p_h_gp = p_h_gp + temp2*R{1}*U(sctr + noDofs*(m-1),:);
        end    
    else
        U_sctr = U(sctr,:);
        p_h_gp = R{1}*U_sctr;
    end

    if exteriorSHBC
        dp_h_gp = -dp_inc(Y,normals); 
    else
        if d_p == 2
            dp_h_gp = rho*omega^2*( normals(:,1).*(R{1}*Ux(elementSolid(e,:),:)) ...
                                   +normals(:,2).*(R{1}*Uy(elementSolid(e,:),:)) ...
                                   +normals(:,3).*(R{1}*Uz(elementSolid(e,:),:)));
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



