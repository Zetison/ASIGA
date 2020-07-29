function F = applyHWBC(varCol,no_angles)
dp_inc = varCol.dp_inc;
noDofs_tot = varCol.noDofs_tot;
weights = varCol.weights;
controlPts = varCol.controlPts;
degree = varCol.degree(1:2);
elRange = varCol.elRange;
d_p = varCol.patches{1}.nurbs.d_p;
knotVecs = varCol.knotVecs;

[zeta0Nodes, noElems, element, element2, index, pIndex, n_en] = meshBoundary(varCol,0);

Fvalues = zeros(n_en,noElems,no_angles);
indices = zeros(n_en,noElems);

[Q, W] = gaussTensorQuad(degree+1);

parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(d_p-1,2);
    for i = 1:d_p-1
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctrXiEta = zeta0Nodes(element(e,:));
    pts = controlPts(sctrXiEta,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:); % New
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^(d_p-1);
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,d_p-1);
    X = R{1}*pts;
    n = -crossProd./repmat(J_1,1,3);
    deriv = -dp_inc(X,n);
    indices(:,e) = sctrXiEta';    
    Fvalues(:,e,:) = R{1}'*(deriv.*J_1 * J_2 .* W);
end

F = zeros(noDofs_tot,no_angles);        % external force vector
for i = 1:no_angles
    F(:,i) = vectorAssembly(Fvalues(:,:,i),indices,noDofs_tot);
end
