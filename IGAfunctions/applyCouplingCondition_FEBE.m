function A = applyCouplingCondition_FEBE(varCol)

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
noDofs = varCol.noCtrlPts;

d = 3;

n_en = (p_xi+1)*(p_eta+1);
A1values = zeros(d*n_en^2,noElems);

spIdxRow1 = zeros(d*n_en^2,noElems);
spIdxCol1 = zeros(d*n_en^2,noElems);


[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
Qxi = Q2D(:,1);
Qeta = Q2D(:,2);
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:); % New

    xi   = parent2ParametricSpace(Xi_e,  Qxi);
    eta  = parent2ParametricSpace(Eta_e, Qeta);
    [R, dRdxi, dRdeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
    J1 = dRdxi*pts;
    J2 = dRdeta*pts;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    normal = crossProd./J_1;

    A1_e = zeros(n_en,d*n_en);
    for i = 1:numel(W2D)
        A1_e = A1_e + R(i,:)'*kron(R(i,:),normal(i,:))*norm(crossProd(i,:)) * J_2 * W2D(i); 
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
