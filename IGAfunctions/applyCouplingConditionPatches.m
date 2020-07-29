function A = applyCouplingConditionPatches(varColInner,varColOuter,shift,shift2,noDofs_tot)
% Future work: Vectorize over Gauss points

knotVecs = varColInner.knotVecs;
elRange = varColInner.elRange;
d_p = varColInner.patches{1}.nurbs.d_p;

degree = varColInner.degree(1:2);


weights = varColInner.weights;
controlPts = varColInner.controlPts;

d_inner = varColInner.dimension;
d_outer = varColOuter.dimension;

[innerNodes,outerNodes,innerXiEtaMesh,innerIndexXiEta,innerNoElemsXiEta,pIndex,innerNodes2] ...
                    = createSurfaceMesh(varColInner,varColOuter);

n_en = prod(degree+1);
Avalues = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);

spIdxRow1 = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);
spIdxCol1 = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);


[Q, W] = gaussTensorQuad(degree+1);
parfor e = 1:innerNoElemsXiEta
% for e = 1:innerNoElemsXiEta    
    patch = pIndex(e); % New
    knots = knotVecs{patch}(1:2);

    Xi_e = zeros(d_p-1,2);
    for i = 1:d_p-1
        Xi_e(i,:) = elRange{i}(innerIndexXiEta(e,i),:);
    end

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^(d_p-1);

    innerSctrXiEta = innerNodes(innerXiEtaMesh(e,:));          %  element scatter vector
    innerSctrXiEtadD = zeros(d_inner*length(innerSctrXiEta),1);
    for i = 1:d_inner
        innerSctrXiEtadD(i:d_inner:d_inner*n_en) = d_inner*(innerSctrXiEta-1)+i;
    end
    outerSctrXiEta = outerNodes(innerXiEtaMesh(e,:));          %  element scatter vector
    outerSctrXiEtadD = zeros(d_outer*length(outerSctrXiEta),1);
    for i = 1:d_outer
        outerSctrXiEtadD(i:d_outer:d_outer*n_en) = d_outer*(outerSctrXiEta-1)+i; % :d_outer
    end

    pts = controlPts(innerSctrXiEta,:);
    wgts = weights(innerNodes2(innerXiEtaMesh(e,:)));
        
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,2);
    normal = crossProd./repmat(J_1,1,3);
    spIdxRow1(:,e) = copyVector(outerSctrXiEtadD,d_inner*n_en,1);
    spIdxCol1(:,e) = copyVector(innerSctrXiEtadD,d_outer*n_en,2);
    A_e = zeros(d_outer*n_en,d_inner*n_en);
    for gp = 1:size(W,1)
        if d_outer > d_inner
            A_e = A_e + kron(kron(R{1}(gp,:),normal(gp,:)),R{1}(gp,:)').'*J_1(gp) * J_2 * W(gp);       
        else
            A_e = A_e + kron(kron(R{1}(gp,:),normal(gp,:)),R{1}(gp,:)')*J_1(gp) * J_2 * W(gp);       
        end
    end
    Avalues(:,e) = reshape(A_e,d_inner*d_outer*n_en^2,1);   
    
    spIdxRow1(:,e) = copyVector(outerSctrXiEtadD,d_inner*n_en,1);
    spIdxCol1(:,e) = copyVector(innerSctrXiEtadD,d_outer*n_en,2);
end

spIdxRow1 = reshape(spIdxRow1,numel(spIdxRow1),1);
spIdxCol1 = reshape(spIdxCol1,numel(spIdxCol1),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx1,~,IuniqueIdx1] = unique([spIdxRow1+shift, spIdxCol1+shift2],'rows');
Avalues = accumarray(IuniqueIdx1,Avalues);

A = sparse(spIdx1(:,1),spIdx1(:,2),Avalues,noDofs_tot,noDofs_tot,numel(IuniqueIdx1));
A = A + A.';
