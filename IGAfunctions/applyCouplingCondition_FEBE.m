function A = applyCouplingCondition_FEBE(varCol)

degree = varCol.degree(1:2);
elRange = varCol.elRange;

index = varCol.index;
noElems = varCol.noElems;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
noDofs = varCol.noCtrlPts;

d = 3;

n_en = prod(degree+1);
A1values = zeros(d*n_en^2,noElems);

spIdxRow1 = zeros(d*n_en^2,noElems);
spIdxCol1 = zeros(d*n_en^2,noElems);

[Q, W] = gaussTensorQuad(degree+1);
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e); % New
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(2,2);
    for i = 1:2
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;

    sctr = element(e,:);
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:); % New
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,2);
    normal = crossProd./J_1;

    A1_e = zeros(n_en,d*n_en);
    for i = 1:numel(W)
        A1_e = A1_e + R{1}(i,:)'*kron(R{1}(i,:),normal(i,:))*norm(crossProd(i,:)) * J_2 * W(i); 
    end
    
    spIdxRow1(:,e) = copyVector(sctr,d*n_en,1);
    spIdxCol1(:,e) = copyVector(sctr_k_e,n_en,2);
    
    A1values(:,e) = reshape(A1_e,d*n_en^2,1);
end

spIdxRow1 = reshape(spIdxRow1,numel(spIdxRow1),1);
spIdxCol1 = reshape(spIdxCol1,numel(spIdxCol1),1);
A1values = reshape(A1values,numel(A1values),1);

[spIdx1,~,IuniqueIdx1] = unique([spIdxRow1, spIdxCol1],'rows');
A1values = accumarray(IuniqueIdx1,A1values);

A = sparse(spIdx1(:,1),spIdx1(:,2),A1values,noDofs,d*noDofs,numel(IuniqueIdx1));
