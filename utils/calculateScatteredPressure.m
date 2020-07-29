function p_h = calculateScatteredPressure(varCol, U, P_far, useExtraQuadPts, computeFarField)

degree = varCol.degree; % assume p_xi is equal in all patches

noDofs = varCol.noDofs;
index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
BC = varCol.BC;
k = varCol.k;
method = varCol.method;
if strcmp(method,'IENSG')
    A_2 = varCol.A_2;
    x_0 = varCol.x_0;
    Upsilon = varCol.Upsilon;
    D = varCol.D;
    N = varCol.N;
else
    A_2 = NaN;
    x_0 = NaN;
    Upsilon = NaN;
    D = NaN;
    N = NaN;
end

dp_inc = varCol.dp_inc;

Phi_k = varCol.Phi_k;
d_p = varCol.patches{1}.nurbs.d_p;

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
solveForPtot = varCol.solveForPtot;
SHBC = strcmp(varCol.BC, 'SHBC');
if SHBC
    if solveForPtot
        homNeumanCond = true;
        dpdn = @(x,n) 0;
    else
        homNeumanCond = false;
        dpdn = @(x,n) -varCol.dp_inc(x,n);
    end
else
    homNeumanCond = false;
    dpdn = varCol.dpdn;
end

extraGP = varCol.extraGP;
if useExtraQuadPts
    noGpXi = degree(1)+1+5;
    noGpEta = degree(2)+1+5;
else
    noGpXi = degree(1)+1+extraGP;
    noGpEta = degree(2)+1+extraGP;
end
[Q, W] = gaussTensorQuad([noGpXi,noGpEta]); 

p_h = zeros(size(P_far,1),1);

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

    if strcmp(BC, 'SHBC')
        dp_h_gp = -dp_inc(Y,normals); 
    else
        if d_p == 2
            error('Not implemented')
        end
        dXdxi = R{2}*pts;
        dXdeta = R{3}*pts;
        dXdzeta = R{4}*pts;
        dp_h_gp = zeros(size(p_h_gp));
        for gp = 1:size(W,1)
            J = [dXdxi(gp,:).' dXdeta(gp,:).' dXdzeta(gp,:).'];
            dp_h_gp(gp,:) = (J'\[R{2}(gp,:); R{3}(gp,:); R{4}(gp,:)]*U_sctr).'*normals(gp,:).';
        end
    end

    X = P_far./repmat(norm2(P_far),1,size(P_far,2));
    
    if computeFarField
        x_d_n = normals*X.';
        x_d_y = Y*X.';
        if solveForPtot
            p_h = p_h + (1i*k*p_h_gp.*x_d_n.*exp(-1i*k*x_d_y)).'* (J_1 * J_2 .* W);  
        else
            p_h = p_h + ((1i*k*p_h_gp.*x_d_n + dp_h_gp).*exp(-1i*k*x_d_y)).'* (J_1 * J_2 .* W);  
        end
    else
        xmy = reshape(P_far,[1,size(P_far,1),3]) - reshape(Y,[size(Y,1),1,3]);
        r = norm2(xmy);
        dPhidny =  Phi_k(r)./r.^2.*(1 - 1i*k*r).*sum(xmy.*repmat(reshape(normals,size(normals,1),1,3),1,size(P_far,1),1),3);
        if solveForPtot
            p_h = p_h + (p_h_gp.*dPhidny).'* (J_1 * J_2 .* W); 
        else 
            p_h = p_h + (p_h_gp.*dPhidny - dp_h_gp.*Phi_k(r)).'* (J_1 * J_2 .* W);
        end 
    end
end
if computeFarField
    p_h = -1/(4*pi)*p_h;
end



