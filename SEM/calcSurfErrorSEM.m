function relError = calcSurfErrorSEM(varCol, LpOrder)
% error('Implementation not complete')
%% Preallocation and initiallizations
patches = varCol.patches;
nurbs = patches{1}.nurbs;
U = varCol{1}.U;

Nxi = nurbs.number(1);
Neta = nurbs.number(2);
Nzeta = nurbs.number(3);
nxi = Nxi+5;
neta = Neta+5;

noComponents = 1;
noPatches = varCol.noPatches;
nzeta = 1;
factors = zeros(nxi*neta*nzeta,noPatches);
nodes = zeros(nxi*neta*nzeta,noPatches,3);
u_hs = zeros(nxi*neta*nzeta,1,noPatches);
[Qxi, Wxi] = gaussTensorQuad(nxi);
[Qeta, Weta] = gaussTensorQuad(neta);

parfor patch = 1:noPatches
% for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    U_sctr = U((1:Nxi*Neta*Nzeta) + (patch-1)*Nxi*Neta*Nzeta,:);

    c = nurbs.coeffs;
    Bxi = zeros(Nxi,nxi);
    Beta = zeros(Neta,neta);
    Bzeta = zeros(Nzeta,1);
    dBxi = zeros(Nxi,nxi);
    dBeta = zeros(Neta,neta);
    GLL = nurbs.GLL;
    for i = 1:Nxi
        Bxi(i,:) = lagrangePolynomials(Qxi,i,Nxi,GLL{1});
        dBxi(i,:) = lagrangePolynomialsDeriv(Qxi,i,Nxi,GLL{1});
    end
    for j = 1:Neta
        Beta(j,:) = lagrangePolynomials(Qeta,j,Neta,GLL{2});
        dBeta(j,:) = lagrangePolynomialsDeriv(Qeta,j,Neta,GLL{2});
    end
    Bzeta(1,:) = 1;
    % nxi = length(xi);
    % neta = length(eta);

    XX = zeros(nxi,neta,1,3);
    JJ = zeros(nxi,neta,1,3,2);
    idx = [2,3,1];
    u_h =           permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(U_sctr,Nxi,Neta*Nzeta)), nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);
    for i = 1:3
        XX(:,:,:,i) = permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);
        JJ(:,:,:,i,1) = permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((dBxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
        JJ(:,:,:,i,2) = permute(reshape(Bzeta.'*reshape(permute(reshape(dBeta.'*reshape(permute(reshape((Bxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
    end
    WW = copyVector2(Wxi.',neta,1,1).*copyVector2(Weta.',nxi,1,2);
    JJ = reshape(permute(JJ,[4,5,1,2,3]),3,2,nxi*neta);
    J = norm2(reshape(cross(JJ(:,1,:),JJ(:,2,:),1),3,[]).');
    
    factors(:,patch) = J.*WW;
    nodes(:,patch,:) = reshape(XX, nxi*neta,1,3);
    u_hs(:,patch) = reshape(u_h, nxi*neta,noComponents);
end
factors = reshape(factors,nxi*neta*noPatches,1);
nodes = reshape(nodes,nxi*neta*noPatches,3);

p_h = reshape(u_hs,nxi*neta*noPatches,1);

analyticFunctions = varCol.analyticFunctions({nodes});
p = analyticFunctions{1}.p;

if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*factors))^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*factors))^(1/LpOrder);
end

relError = 100*Error/normalization;



function d = det2(A)

d = A(1,1,:,:,:).*A(2,2,:,:,:) - A(1,2,:,:,:).*A(2,1,:,:,:);

function d = det3(A)

d = A(1,1,:,:,:).*det2(A(2:3,2:3,:,:,:)) - A(1,2,:,:,:).*det2(A(2:3,[1,3],:,:,:)) + A(1,3,:,:,:).*det2(A(2:3,1:2,:,:,:));

