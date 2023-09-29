function p_h = calculateScatteredPressure(task, P_far, useExtraQuadPts)

plotFarField = task.ffp.plotFarField;
d_p = task.varCol{1}.patches{1}.nurbs.d_p;
d = task.varCol{1}.patches{1}.nurbs.d;
noDofs = task.varCol{1}.noDofs;
noElems = task.varCol{1}.noElems;
degree = task.varCol{1}.degree(1:2); % assume p_xi is equal in all patches

weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;
rho = task.varCol{1}.rho;
if d_p == 3
    varColBdry = meshBoundary(task.varCol{1},'Gamma');
    
    zeta0Nodes = varColBdry.nodes;
    noElems = varColBdry.noElems;
    element = varColBdry.element;
    element2 = varColBdry.element2;
    index = varColBdry.index;
    pIndex = varColBdry.pIndex;
    knotVecs = varColBdry.knotVecs;
    elRange = varColBdry.elRange;
    n_en = varColBdry.n_en;
else
    degree = task.varCol{1}.degree; % assume p_xi is equal in all patches
    index = task.varCol{1}.index;
    element = task.varCol{1}.element;
    element2 = task.varCol{1}.element2;
    pIndex = task.varCol{1}.pIndex;
    zeta0Nodes = 1:noDofs;
    knotVecs = task.varCol{1}.knotVecs;
    elRange = task.varCol{1}.elRange;
    n_en = prod(degree+1);
end

BC = task.misc.BC;
omega = task.misc.omega;
k = omega/task.varCol{1}.c_f;

dp_inc = task.dp_incdn_;

Phi_k = task.varCol{1}.Phi_k;
solveForPtot = task.misc.solveForPtot;
exteriorSHBC = (strcmp(BC, 'SHBC') || strcmp(BC, 'NBC')) && numel(task.varCol) == 1;

extraGP = max(task.misc.extraGP,task.ffp.extraGP);
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

noRHSs = size(task.varCol{1}.U,2);
Uvalues = zeros(n_en,noElems,noRHSs);
for e = 1:noElems
    sctr = zeta0Nodes(element(e,:));
    Uvalues(:,e,:) = task.varCol{1}.U(sctr,:);
end
if numel(task.varCol) > 1
    varColBdrySolid = meshBoundary(task.varCol{2},'Gamma');
    
    nodesSolid = varColBdrySolid.nodes;
    elementSolid = varColBdrySolid.element;   

    Uxvalues = zeros(n_en,noElems,noRHSs);
    Uyvalues = zeros(n_en,noElems,noRHSs);
    Uzvalues = zeros(n_en,noElems,noRHSs); 

    for e = 1:noElems
        sctrSolid = nodesSolid(elementSolid(e,:))*d;
        Uxvalues(:,e,:) = task.varCol{2}.U(sctrSolid-2,:);
        Uyvalues(:,e,:) = task.varCol{2}.U(sctrSolid-1,:);
        Uzvalues(:,e,:) = task.varCol{2}.U(sctrSolid,:);
    end
else
    Uxvalues = NaN(1,noElems);
    Uyvalues = NaN(1,noElems);
    Uzvalues = NaN(1,noElems);
end
% for e = 1:noElems
parfor e = 1:noElems
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
    p_h_gp = R{1}*squeeze(Uvalues(:,e,:));

    if exteriorSHBC
        dp_h_gp = -dp_inc(Y,normals); 
    else
        dp_h_gp = rho*omega.^2.*( normals(:,1).*(R{1}*squeeze(Uxvalues(:,e,:))) ...
                                 +normals(:,2).*(R{1}*squeeze(Uyvalues(:,e,:))) ...
                                 +normals(:,3).*(R{1}*squeeze(Uzvalues(:,e,:))));
        if ~solveForPtot
            dp_h_gp = dp_h_gp - dp_inc(Y,normals);
        end
    end

    X = P_far./repmat(norm2(P_far),1,size(P_far,2));
    
    if plotFarField
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
if plotFarField
    p_h = -1/(4*pi)*p_h;
end



