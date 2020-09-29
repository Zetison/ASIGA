function F = applyHWBC_ROM(varCol,noVecs)

p_inc = varCol.p_inc;
d_vec = varCol.d_vec;
k = varCol.k;
noDofs_tot = varCol.noDofs_tot;
weights = varCol.weights;
controlPts = varCol.controlPts;
degree = varCol.degree(1:2);
elRange = varCol.elRange;
d_p = varCol.patches{1}.nurbs.d_p;
knotVecs = varCol.knotVecs;

[zeta0Nodes, noElems, element, element2, index, pIndex, n_en] = meshBoundary(varCol,0);

Fvalues = zeros(n_en,noElems,noVecs);
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
    
    f_e = zeros(n_en,noVecs);
    for gp = 1:size(W,1)
        m = 0:(noVecs-1);
        d_vecX = dot(d_vec,X(gp,:));
        temp = m/d_vecX;
        temp(1) = 0;
        deriv = -exp(1i*pi*m/2)*dot(d_vec,n(gp,:))*p_inc(X(gp,:)).*d_vecX.^m.*(temp+1i*k);
        f_e = f_e + R{1}(gp,:)'*deriv*J_1(gp) * J_2 * W(gp);  
    end  
    indices(:,e) = sctrXiEta';
    Fvalues(:,e,:) = f_e;
end

F = zeros(noDofs_tot,noVecs);        % external force vector
for i = 1:noVecs
    F(:,i) = vectorAssembly(Fvalues(:,:,i),indices,noDofs_tot);
end

