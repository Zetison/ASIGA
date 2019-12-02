function F = applyNeumannCondition_CoupledProblemPatches(varColSolid,varColFluid,omega,rho_f,no_angles,shift)

knotVecs = varColSolid.knotVecs;
elRangeXi = varColSolid.elRange{1};
elRangeEta = varColSolid.elRange{2};

p_xi = varColSolid.degree(1);
p_eta = varColSolid.degree(2);

weights = varColSolid.weights;
controlPts = varColSolid.controlPts;

dp_inc = varColSolid.dp_inc;
p_inc = varColSolid.p_inc;


noDofs = varColSolid.noDofs;
noDofs_tot = varColSolid.noDofs_tot;
d = varColSolid.dimension;

[solidNodes,fluidNodes,solidXiEtaMesh,solidIndexXiEta,solidNoElemsXiEta,pIndex,solidNodes2,fluidNodes2]...
    = createSurfaceMesh(varColSolid,varColFluid);

n_en = (p_xi+1)*(p_eta+1);

F1values = zeros(d*n_en,solidNoElemsXiEta,no_angles);
F2values = zeros(n_en,solidNoElemsXiEta,no_angles);

indices1 = zeros(d*n_en,solidNoElemsXiEta);
indices2 = zeros(n_en,solidNoElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
% warning('parfor is not being used')
for e = 1:solidNoElemsXiEta
% parfor e = 1:solidNoElemsXiEta
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    
    idXi = solidIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = solidIndexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
    solidSctrXiEtadD = zeros(d*length(solidSctrXiEta),1);
    for i = 1:d
        solidSctrXiEtadD(i:d:d*n_en) = d*(solidSctrXiEta-1)+i;
    end
    fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:))+noDofs;          %  element scatter vector

    pts = controlPts(solidSctrXiEta,:);
    wgts = weights(solidNodes2(solidXiEtaMesh(e,:)));
    
    F1_e = zeros(d*n_en,no_angles);
    F2_e = zeros(n_en,no_angles);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        J_1 = norm(crossProd);
        n = crossProd/J_1;
        
        X = R*pts;
                
        F1_e = F1_e + kron(R',n)*p_inc(X).'*J_1 * J_2 * wt;
        F2_e = F2_e + R'*dp_inc(X,n.').'*J_1 * J_2 * wt;
    end
    
    indices1(:,e) = solidSctrXiEtadD';
    indices2(:,e) = fluidSctrXiEta';
    F1values(:,e,:) = F1_e;
    F2values(:,e,:) = F2_e;
end

for alpha_s_Nr = 1:no_angles
    F(:,alpha_s_Nr) = -vectorAssembly(F1values(:,:,alpha_s_Nr),indices1+shift,noDofs_tot);
    F(:,alpha_s_Nr) = F(:,alpha_s_Nr) + 1/(rho_f*omega^2)*vectorAssembly(F2values(:,:,alpha_s_Nr),indices2+shift,noDofs_tot);
end

