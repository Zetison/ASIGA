function [relL2Error, relH1Error, relH1sError, relEnergyError] = calcErrorSEM(varCol, UU, options)

%% Preallocation and initiallizations
patches = varCol.patches;
nurbs = patches{1}.nurbs;

Nxi = nurbs.number(1);
Neta = nurbs.number(2);
Nzeta = nurbs.number(3);
nxi = Nxi+5;
neta = Neta+5;
nzeta = Nzeta+5;

noComponents = 1;
noComponentsDeriv = 3;
noPatches = varCol.noPatches;
factors = zeros(nxi*neta*nzeta,noPatches);
nodes = zeros(nxi*neta*nzeta,noPatches,3);
u_hs = zeros(nxi*neta*nzeta,1,noPatches);
du_hs = zeros(nxi*neta*nzeta,noPatches,3);
[Wxi,Qxi] = gaussianQuadNURBS(nxi); 
[Weta,Qeta] = gaussianQuadNURBS(neta); 
[Wzeta,Qzeta] = gaussianQuadNURBS(nzeta); 

parfor patch = 1:noPatches
% for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    U = UU((1:Nxi*Neta*Nzeta) + (patch-1)*Nxi*Neta*Nzeta,:);

    c = nurbs.coeffs;
    Bxi = zeros(Nxi,nxi);
    Beta = zeros(Neta,neta);
    Bzeta = zeros(Nzeta,nzeta);
    dBxi = zeros(Nxi,nxi);
    dBeta = zeros(Neta,neta);
    dBzeta = zeros(Nzeta,nzeta);
    GLL = nurbs.GLL;
    for i = 1:Nxi
        Bxi(i,:) = lagrangePolynomials(Qxi,i,Nxi,GLL{1});
        dBxi(i,:) = lagrangePolynomialsDeriv(Qxi,i,Nxi,GLL{1});
    end
    for j = 1:Neta
        Beta(j,:) = lagrangePolynomials(Qeta,j,Neta,GLL{2});
        dBeta(j,:) = lagrangePolynomialsDeriv(Qeta,j,Neta,GLL{2});
    end
    for l = 1:Nzeta
        Bzeta(l,:) = lagrangePolynomials(Qzeta,l,Nzeta,GLL{3});
        dBzeta(l,:) = lagrangePolynomialsDeriv(Qzeta,l,Nzeta,GLL{3});
    end
    % nxi = length(xi);
    % neta = length(eta);

    XX = zeros(nxi,neta,nzeta,3);
    du_ht = zeros(nxi,neta,nzeta,noComponentsDeriv);
    JJ = zeros(nxi,neta,nzeta,3,3);
    idx = [2,3,1];
    u_h =           permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(U,Nxi,Neta*Nzeta)), nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);
    du_ht(:,:,:,1) = permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((dBxi.'*reshape(U,Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
    du_ht(:,:,:,2) = permute(reshape(Bzeta.'*reshape(permute(reshape(dBeta.'*reshape(permute(reshape((Bxi.'*reshape(U,Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
    du_ht(:,:,:,3) = permute(reshape(dBzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(U,Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
    for i = 1:3
        XX(:,:,:,i) = permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);
        JJ(:,:,:,i,1) = permute(reshape(Bzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((dBxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
        JJ(:,:,:,i,2) = permute(reshape(Bzeta.'*reshape(permute(reshape(dBeta.'*reshape(permute(reshape((Bxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
        JJ(:,:,:,i,3) = permute(reshape(dBzeta.'*reshape(permute(reshape(Beta.'*reshape(permute(reshape((Bxi.'*reshape(c(i,:,:,:),Nxi,Neta*Nzeta)),nxi,Neta,Nzeta),idx),Neta,Nzeta*nxi),neta,Nzeta,nxi),idx),Nzeta,nxi*neta),nzeta,nxi,neta),idx);        
    end
    WW = copyVector2(Wxi.',neta,nzeta,1).*copyVector2(Weta.',nxi,nzeta,2).*copyVector2(Wzeta.',nxi,neta,3);
    JJ = permute(JJ,[4,5,1,2,3]);
    J = det3(JJ);
    Jinv = permute(inv3(JJ,J),[3,4,5,1,2]);
    du_h = zeros(nxi,neta,nzeta,noComponentsDeriv);
    du_h(:,:,:,1) = sum(Jinv(:,:,:,:,1).*du_ht,4);
    du_h(:,:,:,2) = sum(Jinv(:,:,:,:,2).*du_ht,4);
    du_h(:,:,:,3) = sum(Jinv(:,:,:,:,3).*du_ht,4);
    
    factors(:,patch) = reshape(J,nxi*neta*nzeta,1).*WW;
    nodes(:,patch,:) = reshape(XX, nxi*neta*nzeta,1,3);
    u_hs(:,patch) = reshape(u_h, nxi*neta*nzeta,noComponents);
    du_hs(:,patch,:) = reshape(du_h, nxi*neta*nzeta,1,noComponentsDeriv);
end
factors = reshape(factors,nxi*neta*nzeta*noPatches,1);
nodes = reshape(nodes,nxi*neta*nzeta*noPatches,3);

u_hs = reshape(u_hs,nxi*neta*nzeta*noPatches,1);
du_hs = reshape(du_hs,nxi*neta*nzeta*noPatches,3);
if strcmp(varCol.applyLoad, 'radialPulsation')
    data.p = varCol.analytic(nodes);
    dp = varCol.gAnalytic(nodes);
    data.dpdx = dp(:,1);
    data.dpdy = dp(:,2);
    data.dpdz = dp(:,3);
else
    data = e3Dss(nodes,options);
end
rho_f = options.rho_f;
omega = options.omega;
c_f = options.c_f;
k = omega./c_f;

H1Error = 0;
H1sError = 0;
EnergyError = 0;
L2Error = 0;
normalizationEnergy = 0;
normalizationH1 = 0;
normalizationH1s = 0;
normalizationL2 = 0;
m = 1;
p = data(m).p;
dp = [data(m).dpdx, data(m).dpdy, data(m).dpdz];
p_e = p-u_hs;
dp_e = dp-du_hs;

p2 = abs(p).^2;
dp2 = sum(abs(dp).^2,2);
p_e2 = abs(p_e).^2;
dp_e2 = sum(abs(dp_e).^2,2);

H1Error = H1Error + sum((p_e2 + dp_e2).*factors);
H1sError = H1sError + sum(dp_e2.*factors);
normalizationH1 = normalizationH1 + sum((p2 + dp2).*factors);
normalizationH1s = normalizationH1s + sum(dp2.*factors);

EnergyError = EnergyError + 1/(rho_f(m)*omega^2)*sum((dp_e2 + k^2*p_e2).*factors);
normalizationEnergy = normalizationEnergy + 1/(rho_f(m)*omega^2)*sum((dp2 + k^2*p2).*factors);

L2Error = L2Error + sum(p_e2.*factors);
normalizationL2 = normalizationL2 + sum(p2.*factors);

relEnergyError = 100*sqrt(EnergyError/normalizationEnergy);
relH1Error = 100*sqrt(H1Error/normalizationH1);
relH1sError = 100*sqrt(H1sError/normalizationH1s);
relL2Error = 100*sqrt(L2Error/normalizationL2);


function d = det2(A)

d = A(1,1,:,:,:).*A(2,2,:,:,:) - A(1,2,:,:,:).*A(2,1,:,:,:);

function d = det3(A)

d = A(1,1,:,:,:).*det2(A(2:3,2:3,:,:,:)) - A(1,2,:,:,:).*det2(A(2:3,[1,3],:,:,:)) + A(1,3,:,:,:).*det2(A(2:3,1:2,:,:,:));

