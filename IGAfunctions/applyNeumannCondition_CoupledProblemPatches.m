function varColFluid = applyNeumannCondition_CoupledProblemPatches(varColSolid,varColFluid,omega,rho_f,no_angles,shift)

knotVecs = varColSolid.knotVecs;
degree = varColSolid.degree(1:2);
elRange = varColSolid.elRange;
d_p = varColSolid.patches{1}.nurbs.d_p;

weights = varColSolid.weights;
controlPts = varColSolid.controlPts;

dp_inc = varColSolid.dp_inc;
p_inc = varColSolid.p_inc;


noDofs = varColSolid.noDofs;
noDofs_tot = varColSolid.noDofs_tot;
d_f = varColSolid.fieldDimension;

[solidNodes,fluidNodes,solidXiEtaMesh,solidIndexXiEta,solidNoElemsXiEta,pIndex,solidNodes2]...
    = createSurfaceMesh(varColSolid,varColFluid);

n_en = prod(degree+1);

F1values = zeros(d_f*n_en,solidNoElemsXiEta,no_angles);
F2values = zeros(n_en,solidNoElemsXiEta,no_angles);

indices1 = zeros(d_f*n_en,solidNoElemsXiEta);
indices2 = zeros(n_en,solidNoElemsXiEta);

[Q, W] = gaussTensorQuad(degree+1);
% warning('parfor is not being used')
% for e = 1:solidNoElemsXiEta
parfor e = 1:solidNoElemsXiEta
    patch = pIndex(e); % New
    knots = knotVecs{patch}(1:2);

    Xi_e = zeros(d_p-1,2);
    for i = 1:d_p-1
        Xi_e(i,:) = elRange{i}(solidIndexXiEta(e,i),:);
    end

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^(d_p-1);

    solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
    solidSctrXiEtadD = zeros(d_f*length(solidSctrXiEta),1);
    for i = 1:d_f
        solidSctrXiEtadD(i:d_f:d_f*n_en) = d_f*(solidSctrXiEta-1)+i;
    end
    fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:))+noDofs;          %  element scatter vector

    pts = controlPts(solidSctrXiEta,:);
    wgts = weights(solidNodes2(solidXiEtaMesh(e,:)));
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,2);
    normal = crossProd./repmat(J_1,1,3);

    X = R{1}*pts;
    
    F1_e = zeros(d_f*n_en,no_angles);
    F2_e = zeros(n_en,no_angles);
    for gp = 1:size(W,1)
        F1_e = F1_e + kron(R{1}(gp,:)',normal(gp,:).')*p_inc(X(gp,:)).'*J_1(gp) * J_2 * W(gp);
        F2_e = F2_e + R{1}(gp,:)'*dp_inc(X(gp,:),normal(gp,:)).'*J_1(gp) * J_2 * W(gp);
    end
    
    indices1(:,e) = solidSctrXiEtadD';
    indices2(:,e) = fluidSctrXiEta';
    F1values(:,e,:) = F1_e;
    F2values(:,e,:) = F2_e;
end
F = zeros(noDofs_tot,no_angles);
for alpha_s_Nr = 1:no_angles
    F(:,alpha_s_Nr) = -vectorAssembly(F1values(:,:,alpha_s_Nr),indices1+shift,noDofs_tot);
    F(:,alpha_s_Nr) = F(:,alpha_s_Nr) + 1/(rho_f*omega^2)*vectorAssembly(F2values(:,:,alpha_s_Nr),indices2+shift,noDofs_tot);
end
varColFluid.FF = F;
