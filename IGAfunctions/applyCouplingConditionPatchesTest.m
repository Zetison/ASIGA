function A = applyCouplingConditionPatchesTest(varColInner,varColOuter,shift,shift2,noDofs_tot)

knotVecs = varColInner.knotVecs;
elRangeXi = varColInner.elRange{1};
elRangeEta = varColInner.elRange{2};

p_xi = varColInner.degree(1);
p_eta = varColInner.degree(2);


weights = varColInner.weights;
controlPts = varColInner.controlPts;

noDofs = varColInner.noDofs;

d_inner = varColInner.dimension;
d_outer = varColOuter.dimension;

[innerNodes,outerNodes,innerXiEtaMesh,innerIndexXiEta,innerNoElemsXiEta,pIndex,innerNodes2] ...
                    = createSurfaceMesh(varColInner,varColOuter);

n_en = (p_xi+1)*(p_eta+1);
Avalues = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);

spIdxRow1 = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);
spIdxCol1 = zeros(d_inner*d_outer*n_en^2,innerNoElemsXiEta);


[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
parfor e = 1:innerNoElemsXiEta
% for e = 1:innerNoElemsXiEta
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    
    idXi = innerIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = innerIndexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    innerSctrXiEta = innerNodes(innerXiEtaMesh(e,:));          %  element scatter vector
    innerSctrXiEtadD = zeros(d_inner*length(innerSctrXiEta),1);
    for i = 1:d_inner
        innerSctrXiEtadD(i:d_inner:d_inner*n_en) = d_inner*(innerSctrXiEta-1)+i;
    end
    outerSctrXiEta = outerNodes(innerXiEtaMesh(e,:))+shift2;          %  element scatter vector
    outerSctrXiEtadD = zeros(d_outer*length(outerSctrXiEta),1);
    for i = 1:d_outer
        outerSctrXiEtadD(i:d_outer:d_outer*n_en) = d_outer*(outerSctrXiEta-1)+i;
    end

    pts = controlPts(innerSctrXiEta,:);
    wgts = weights(innerNodes2(innerXiEtaMesh(e,:)));
    
    A_e = zeros(n_en,d_inner*n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        
        normal = crossProd/norm(crossProd);
        
        A_e = A_e + kron(kron(R,normal'),R')*norm(crossProd) * J_2 * wt;      
    end
    spIdxRow1(:,e) = copyVector(outerSctrXiEtadD,d_inner*n_en,1);
    spIdxCol1(:,e) = copyVector(innerSctrXiEtadD,d_outer*n_en,2);
    
%     temp = zeros(n_en,d_inner*n_en);
%     for j = 1:d_inner
%         temp(:, j:d_inner:end) = A_e(:, 1+(j-1)*n_en:j*n_en);
%     end
    Avalues(:,e) = reshape(A_e,d_inner*d_outer*n_en^2,1);    
end

spIdxRow1 = reshape(spIdxRow1,numel(spIdxRow1),1);
spIdxCol1 = reshape(spIdxCol1,numel(spIdxCol1),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx1,~,IuniqueIdx1] = unique([spIdxRow1, spIdxCol1]+shift,'rows');
Avalues = accumarray(IuniqueIdx1,Avalues);


% Future work: Vectorize over Gauss points
A = sparse(spIdx1(:,1),spIdx1(:,2),Avalues,noDofs_tot,noDofs_tot,numel(IuniqueIdx1));
A = A + A.';
