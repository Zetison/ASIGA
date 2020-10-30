function [A_fluid_o, FF, dofsToRemove] = buildSEMMatrices(varCol)
patches = varCol.patches;
noPatches = varCol.noPatches;
Nxi = patches{1}.nurbs.number(1);
Neta = patches{1}.nurbs.number(2);
Nzeta = patches{1}.nurbs.number(3);
noDofsPatch = Nxi*Neta*Nzeta;
element = varCol.element;
% noDofs = varCol.noDofs;
N = varCol.N;
% formulation = varCol.formulation;
Upsilon = varCol.Upsilon;
D = varCol.D;
Dt = varCol.Dt;
A_2 = varCol.A_2;
x_0 = varCol.x_0;
r_a = varCol.r_a;
noDofs = varCol.noDofs;
noSurfDofs = noPatches*Nxi*Neta;
no_angles = length(varCol.alpha_s);
dp_dinc = varCol.dp_inc;


Dxi = derivativeMatrix(patches{1}.nurbs.GLL{1},Nxi); % assume Nxi to be the same for all patches
if Neta == Nxi
    Deta = Dxi;
else
    Deta = derivativeMatrix(patches{1}.nurbs.GLL{2},Neta); % assume Neta to be the same for all patches
end
if Nzeta == Nxi
    Dzeta = Dxi;
elseif Nzeta == Neta
    Dzeta = Deta;
else
    Dzeta = derivativeMatrix(patches{1}.nurbs.GLL{3},Nzeta); % assume Nzeta to be the same for all patches
end
[BB,elementGamma,elementInf,zeta1Nodes,zeta0Nodes] = addInfElements3_SEM(varCol);

dofsInInfElements = noDofsPatch*N/Nzeta;

spIdxRow = zeros(noDofsPatch^2,noPatches,'uint32');
spIdxCol = zeros(noDofsPatch^2,noPatches,'uint32');
spIdxRowInf = zeros(dofsInInfElements^2,noPatches,'uint32');
spIdxColInf = zeros(dofsInInfElements^2,noPatches,'uint32');
spIdxM = zeros(noDofsPatch,noPatches,'uint32');
Kvalues = zeros(noDofsPatch^2,noPatches); % assume noDofs to be the same for all patches
Fvalues = zeros(Nxi*Neta,noPatches,no_angles); % assume noDofs to be the same for all patches
Findices = zeros(Nxi*Neta,noPatches); % assume noDofs to be the same for all patches
Kinfvalues = zeros(dofsInInfElements^2,noPatches); % assume noDofs to be the same for all patches
Mvalues = zeros(noDofsPatch,noPatches); % assume noDofs to be the same for all patches
for patch = 1:noPatches
% parfor patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    
    Nxi = nurbs.number(1);
    Neta = nurbs.number(2);
    Nzeta = nurbs.number(3);
    rho = nurbs.rho;
    c = nurbs.coeffs(1:3,:,:,:);


    idx = [2 1 3 4];
    dCdxi = permute(reshape(Dxi.'*reshape(permute(c,idx),Nxi,3*Neta*Nzeta),Nxi,3,Neta,Nzeta),idx);
    idx = [3 2 1 4];
    dCdeta = permute(reshape(Deta.'*reshape(permute(c,idx),Neta,3*Nxi*Nzeta),Neta,Nxi,3,Nzeta),idx);
    idx = [4 2 3 1];
    dCdzeta = permute(reshape(Dzeta.'*reshape(permute(c,idx),Nzeta,3*Nxi*Neta),Nzeta,Nxi,Neta,3),idx);

    JJ = zeros(3,3,Nxi,Neta,Nzeta);
    JJ(:,1,:,:,:) = reshape(dCdxi,3,1,Nxi,Neta,Nzeta);
    JJ(:,2,:,:,:) = reshape(dCdeta,3,1,Nxi,Neta,Nzeta);
    JJ(:,3,:,:,:) = reshape(dCdzeta,3,1,Nxi,Neta,Nzeta);

    J = det3(JJ);
    
    Jinv = inv3(JJ,J);
    
    G = zeros(size(Jinv));
    G(1,1,:,:,:) = reshape(sum(reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(1,2,:,:,:) = reshape(sum(reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(1,3,:,:,:) = reshape(sum(reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(2,2,:,:,:) = reshape(sum(reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(2,3,:,:,:) = reshape(sum(reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(3,3,:,:,:) = reshape(sum(reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    
    G(2,1,:,:,:) = reshape(sum(reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(3,1,:,:,:) = reshape(sum(reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(1,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    G(3,2,:,:,:) = reshape(sum(reshape(Jinv(3,:,:,:,:),3,Nxi*Neta*Nzeta).*reshape(Jinv(2,:,:,:,:),3,Nxi*Neta*Nzeta),1),1,1,Nxi,Neta,Nzeta);
    
    G = G.*repmat(J,3,3,1,1,1);
    G = multiplyByWeights(G,rho,Nxi,Neta,Nzeta,3);
    
    J_M = multiplyByWeights(J,rho,Nxi,Neta,Nzeta,1);
    noDofsPatch = Nxi*Neta*Nzeta;
    Mvalues(:,patch) = reshape(J_M,noDofsPatch,1);
    
    
    Kloc = zeros(Nxi,Neta,Nzeta,Nxi,Neta,Nzeta);

    idxXi = [4,2,3,1,5,6];
    idxEta = [1,5,3,4,2,6];
    idxZeta = [1,2,6,4,5,3];
    G11 = reshape(G(1,1,:,:,:),Nxi,Neta,Nzeta);
    G12 = reshape(G(1,2,:,:,:),Nxi,Neta,Nzeta);
    G13 = reshape(G(1,3,:,:,:),Nxi,Neta,Nzeta);
    G22 = reshape(G(2,2,:,:,:),Nxi,Neta,Nzeta);
    G23 = reshape(G(2,3,:,:,:),Nxi,Neta,Nzeta);
    G33 = reshape(G(3,3,:,:,:),Nxi,Neta,Nzeta);

    Dxir = reshape(Dxi,Nxi,1,1,Nxi,1,1);
    Detar = reshape(Deta,1,Neta,1,1,Neta,1);
    Dzetar = reshape(Dzeta,1,1,Nzeta,1,1,Nzeta);
    
    temp =   repmat(permute(G12,idxEta),1,Neta,1,Nxi,1,1)...
              .*repmat(permute(Dxir,idxXi),1,Neta,Nzeta,1,Neta,1).*repmat(Detar,Nxi,1,Nzeta,Nxi,1,1) ...
           + repmat(permute(G12,idxXi),Nxi,1,1,1,Neta,1)...
              .*repmat(Dxir,1,Neta,Nzeta,1,Neta,1).*repmat(permute(Detar,idxEta),Nxi,1,Nzeta,Nxi,1,1);
    for l = 1:Nzeta
        Kloc(:,:,l,:,:,l) = temp(:,:,l,:,:);
    end
    temp =   repmat(permute(G13,idxZeta),1,1,Nzeta,Nxi,1,1)...
              .*repmat(permute(Dxir,idxXi),1,Neta,Nzeta,1,1,Nzeta).*repmat(Dzetar,Nxi,Neta,1,Nxi,1,1) ...
           + repmat(permute(G13,idxXi),Nxi,1,1,1,1,Nzeta)...
              .*repmat(Dxir,1,Neta,Nzeta,1,1,Nzeta).*repmat(permute(Dzetar,idxZeta),Nxi,Neta,1,Nxi,1,1);
    for j = 1:Neta
        Kloc(:,j,:,:,j,:) = Kloc(:,j,:,:,j,:) + temp(:,j,:,:,1,:);
    end
    temp =   repmat(permute(G23,idxZeta),1,1,Nzeta,1,Neta,1)...
              .*repmat(permute(Detar,idxEta),Nxi,1,Nzeta,1,1,Nzeta).*repmat(Dzetar,Nxi,Neta,1,Nxi,1,1) ...
           + repmat(permute(G23,idxEta),1,Neta,1,1,1,Nzeta)...
              .*repmat(Detar,Nxi,1,Nzeta,1,1,Nzeta).*repmat(permute(Dzetar,idxZeta),Nxi,Neta,1,1,Neta,1);
    for i = 1:Nxi
        Kloc(i,:,:,i,:,:) = Kloc(i,:,:,i,:,:) + temp(i,:,:,1,:,:);
    end
    
    temp = sum(repmat(permute(G11,[7,2,3,4,5,6,1]),Nxi,1,1,Nxi,1,1,1) ...
                             .*repmat(permute(Dxir,[1,2,3,7,5,6,4]),1,Neta,Nzeta,Nxi,1,1,1) ...
                             .*repmat(permute(Dxir,[7,2,3,1,5,6,4]),Nxi,Neta,Nzeta,1,1,1,1), 7);
    for j = 1:Neta
        for l = 1:Nzeta
            Kloc(:,j,l,:,j,l) = Kloc(:,j,l,:,j,l) + temp(:,j,l,:,1,1);
        end
    end
    
    temp = sum(repmat(permute(G22,[1,7,3,4,5,6,2]),1,Neta,1,1,Neta,1,1) ...
                             .*repmat(permute(Detar,[1,2,3,4,7,6,5]),Nxi,1,Nzeta,1,Neta,1,1) ...
                             .*repmat(permute(Detar,[1,7,3,4,2,6,5]),Nxi,Neta,Nzeta,1,1,1,1), 7);
    for i = 1:Nxi
        for l = 1:Nzeta
            Kloc(i,:,l,i,:,l) = Kloc(i,:,l,i,:,l) + temp(i,:,l,1,:,1);
        end
    end
    
    temp = sum(repmat(permute(G33,[1,2,7,4,5,6,3]),1,1,Nzeta,1,1,Nzeta,1) ...
                             .*repmat(permute(Dzetar,[1,2,3,4,5,7,6]),Nxi,Neta,1,1,1,Nzeta,1) ...
                             .*repmat(permute(Dzetar,[1,2,7,4,5,3,6]),Nxi,Neta,Nzeta,1,1,1,1), 7);
    for i = 1:Nxi
        for j = 1:Neta
            Kloc(i,j,:,i,j,:) = Kloc(i,j,:,i,j,:) + temp(i,j,:,1,1,:);
        end
    end
%     testStiffnessMatrix
    
    %% Assembly infinite element matrix    
    JJS = JJ(:,1:2,:,:,end);
    
    Xt = reshape(A_2*reshape(c(:,:,:,end)-repmat(x_0.',1,Nxi,Neta),3,Nxi*Neta),3,Nxi,Neta);
    xt = reshape(Xt(1,:,:),Nxi*Neta,1);
    yt = reshape(Xt(2,:,:),Nxi*Neta,1);
    zt = reshape(Xt(3,:,:),Nxi*Neta,1);

    [~, theta, ~, d1, d2] = evaluateProlateCoords2(xt,yt,zt,Upsilon);

    
    temp = zeros(2,3,Nxi*Neta);
    
%     DPDX(0,1,:) = xt.*(d1+d2)./(2*d1*d2);
    temp(1,1,:) = xt.*zt./(d1.*d2.*sqrt(r_a^2-zt.^2));
    temp(2,1,:) = -yt./(xt.^2+yt.^2);
%     DPDX(0,2,:) = yt.*(d1+d2)./(2*d1.*d2);
    temp(1,2,:) = yt.*zt./(d1.*d2.*sqrt(r_a^2-zt.^2));
    temp(2,2,:) = xt./(xt.^2+yt.^2);
%     DPDX(0,3,:) = (zt.*(d1+d2)+Upsilon*(d2-d1))./(2*d1.*d2);
    temp(1,3,:) = 1./sqrt(r_a^2-zt.^2).*(zt.^2./(d1.*d2)+Upsilon*zt.*(d2-d1)./(d1.*d2.*(d1+d2))-1);
%     temp(or(isnan(temp),isinf(temp))) = 0;
    
    DPDX = zeros(2,3,Nxi*Neta);
    DPDX(1,1,:) = sum(temp(1,:,:).*repmat(A_2(:,1).',1,1,Nxi*Neta),2);
    DPDX(1,2,:) = sum(temp(1,:,:).*repmat(A_2(:,2).',1,1,Nxi*Neta),2);
    DPDX(1,3,:) = sum(temp(1,:,:).*repmat(A_2(:,3).',1,1,Nxi*Neta),2);
    DPDX(2,1,:) = sum(temp(2,:,:).*repmat(A_2(:,1).',1,1,Nxi*Neta),2);
    DPDX(2,2,:) = sum(temp(2,:,:).*repmat(A_2(:,2).',1,1,Nxi*Neta),2);
    DPDX(2,3,:) = sum(temp(2,:,:).*repmat(A_2(:,3).',1,1,Nxi*Neta),2);
    
    DPDX = reshape(DPDX,2,3,Nxi,Neta);

    JJ3 = zeros(2,2,Nxi,Neta);
    idx = [2 1 3 4];
    JJ3(1,1,:,:) = sum(permute(DPDX(1,:,:,:),idx).*JJS(:,1,:,:),1);
    JJ3(1,2,:,:) = sum(permute(DPDX(1,:,:,:),idx).*JJS(:,2,:,:),1);
    JJ3(2,1,:,:) = sum(permute(DPDX(2,:,:,:),idx).*JJS(:,1,:,:),1);
    JJ3(2,2,:,:) = sum(permute(DPDX(2,:,:,:),idx).*JJS(:,2,:,:),1);
    J3 = det2(JJ3);
    J3inv = inv2(JJ3,J3);
    J4 = permute(J3inv,[2,1,3,4]); % transpose J3inv

    J5 = zeros(2,3,Nxi,Neta);
    J5(1,1,:,:) = J4(1,1,:,:).*J4(1,1,:,:);
    J5(1,2,:,:) = J4(1,2,:,:).*J4(1,1,:,:);
    J5(1,3,:,:) = J4(1,2,:,:).*J4(1,2,:,:);
    J5(2,1,:,:) = J4(2,1,:,:).*J4(2,1,:,:);
    J5(2,2,:,:) = J4(2,2,:,:).*J4(2,1,:,:);
    J5(2,3,:,:) = J4(2,2,:,:).*J4(2,2,:,:);
    J5 = J5.*repmat(J3,2,3,1,1);
    
    J5 = multiplyByWeightsAtBoundary2(J5,rho,Nxi,Neta,2,3);
    
    theta = reshape(theta,1,1,Nxi,Neta);
    sinTheta = sin(theta);
    cosTheta2 = cos(theta).^2;
    J5(1,:,:,:) = J5(1,:,:,:).*repmat(sinTheta,1,3,1,1);
    J5(2,:,:,:) = J5(2,:,:,:)./repmat(sinTheta,1,3,1,1);
    J5(3,:,:,:) = J5(2,:,:,:);
    J5(3,:,:,:) = J5(3,:,:,:).*repmat(cosTheta2,1,3,1,1);
    
    A = zeros(Nxi,Neta,Nxi,Neta,5);
    
    sinThetaJ3 = sinTheta.*J3;
    sinThetaJ3(or(isnan(sinThetaJ3(:)),isinf(sinThetaJ3(:)))) = 0;
    temp1 = multiplyByWeightsAtBoundary2(sinThetaJ3,rho,Nxi,Neta,1,1);
    temp2 = multiplyByWeightsAtBoundary2(cosTheta2.*sinThetaJ3,rho,Nxi,Neta,1,1);
    for i = 1:Nxi
        for j = 1:Neta
            A(i,j,i,j,1) = temp1(1,1,i,j);
            A(i,j,i,j,3) = temp2(1,1,i,j);
        end
    end
    
    Dxir = reshape(Dxi,Nxi,1,Nxi,1);
    Dxir2 = reshape(Dxi.',Nxi,1,Nxi,1);
    Detar = reshape(Deta,1,Neta,1,Neta);
    Detar2 = reshape(Deta.',1,Neta,1,Neta);
    
    idx = find(sin(theta) < 10*eps);
    if ~isempty(idx)
%         idx = idx(1);
        ii = mod(idx-1,Neta)+1;
        jj = floor((idx-1)/Neta)+1;
        J5(:,:,ii,jj) = 0;
        idx = idx(1);
        ii = ii(1);
        jj = jj(1);
    end
    aa = [2,4,5];
    for iii = 1:3
        a = aa(iii);
        A(:,:,:,:,a) =  repmat(permute(J5(iii,2,:,:),[1,4,3,2]),Nxi,1,1,Neta).*repmat(Dxir,1,Neta,1,Neta).*repmat(Detar2,Nxi,1,Nxi,1) ...
                      + repmat(permute(J5(iii,2,:,:),[3,2,1,4]),1,Neta,Nxi,1).*repmat(Dxir2,1,Neta,1,Neta).*repmat(Detar,Nxi,1,Nxi,1);

        temp = sum(repmat(permute(J5(iii,1,:,:),[1,4,5,2,3]),Nxi,1,1,1) ...
                                 .*repmat(permute(Dxir,[1,2,5,4,3]),1,Neta,1,1) ...
                                 .*repmat(permute(Dxir,[5,2,1,4,3]),1,Neta,1,1), 5);
        for j = 1:Neta
            A(:,j,:,j,a) = A(:,j,:,j,a) + temp(:,j,:,1);
        end    

        temp = sum(repmat(permute(J5(iii,3,:,:),[3,2,1,5,4]),1,Neta,1,1) ...
                                 .*repmat(permute(Detar,[1,2,3,5,4]),Nxi,1,1,1) ...
                                 .*repmat(permute(Detar,[1,5,3,2,4]),Nxi,1,1,1), 5);
        for i = 1:Nxi
            A(i,:,i,:,a) = A(i,:,i,:,a) + temp(i,:,1,:);
        end
    end
    thetaIdx = theta(idx);
    noSp = numel(thetaIdx); % number of singular points
    if ~isempty(idx)
        J = JJS(:,:,ii,jj);
        sgn = sign(cos(thetaIdx(1)));
        detrm = J(1,1)*J(2,2)-J(1,2)*J(2,1);
        wt = rho{1}(ii(1)).*rho{2}(jj(1));
        A(ii,jj,ii,jj,[1,3]) = A(ii,jj,ii,jj,[1,3]) + sgn*detrm/(r_a^2-Upsilon^2)*wt;
        
        A(:,jj,ii,:,2) = A(:,jj,ii,:,2) - sgn*(J(1,1)*J(1,2)+J(2,1)*J(2,2))*reshape(Dxi(:,ii)*Deta(:,jj).',Nxi,1,1,Neta)/detrm*wt;
        A(ii,:,:,jj,2) = A(ii,:,:,jj,2) - sgn*(J(1,1)*J(1,2)+J(2,1)*J(2,2))*reshape(Deta(:,jj)*Dxi(:,ii).',1,Neta,Nxi,1)/detrm*wt;

        temp = sgn*(J(1,2)^2+J(2,2)^2)*Dxi(:,ii)*Dxi(:,ii).'/detrm*wt;
        A(:,jj,:,jj,2) = A(:,jj,:,jj,2) + reshape(temp,Nxi,1,Nxi,1);

        temp = sgn*(J(1,1)^2+J(2,1)^2)*Deta(:,jj)*Deta(:,jj).'/detrm*wt;
        A(ii,:,ii,:,2) = A(ii,:,ii,:,2) + reshape(temp,1,Neta,1,Neta);
    end
    
    Kinf = zeros(Nxi,Neta,Nxi,Neta,N,N);
    for m = 1:N
        for n = 1:N
            for mt = 1:N 
                for nt = 1:N
                    temp = zeros(Nxi,Neta,Nxi,Neta);
                    for a = 1:5
                        temp = temp + A(:,:,:,:,a)*BB(nt,mt,a);
                    end
                    Kinf(:,:,:,:,n,m) = Kinf(:,:,:,:,n,m) + temp*D(m,mt)*Dt(n,nt); 
                end
            end
        end
    end
    
%     testInfMatrix
    
    sctrInf = elementInf(patch,:);
    sctrGamma = zeta0Nodes(elementGamma(patch,:));
    sctrGlobal = zeta1Nodes(sctrInf);
    spIdxRow_temp = zeros((N*Nxi*Neta)^2,1);
    spIdxCol_temp = zeros((N*Nxi*Neta)^2,1);
    n_en = Nxi*Neta;
    for m = 1:N
        for n = 1:N
            indices = (1:n_en^2)+(n_en^2*N*(m-1) + n_en^2*(n-1));
            if n == 1
                IEsctr = sctrGlobal;
            else
                IEsctr = noDofs+sctrInf+noSurfDofs*(n-2);
            end
            spIdxRow_temp(indices) = copyVector(IEsctr,n_en,1);
            if m == 1
                IEsctr = sctrGlobal;
            else
                IEsctr = noDofs+sctrInf+noSurfDofs*(m-2);
            end
            spIdxCol_temp(indices) = copyVector(IEsctr,n_en,2);
        end
    end
    spIdxRowInf(:,patch) = spIdxRow_temp;
    spIdxColInf(:,patch) = spIdxCol_temp;
    Kinfvalues(:,patch) = reshape(Kinf,dofsInInfElements^2,1);
    
    %% Build RHS
    normals = reshape(cross(dCdxi(:,:,:,1),dCdeta(:,:,:,1),1),3,Nxi*Neta);
    H = reshape(norm2(normals.'),Nxi,Neta);
    
    g = reshape(-dp_dinc(reshape(c(:,:,:,1),3,Nxi*Neta).', -(normals./repmat(reshape(H,1,Nxi*Neta),3,1)).'),Nxi,Neta);
    F = multiplyByWeightsAtBoundary(g.*H,rho,Nxi,Neta);
    
    %% Collect matrices
    sctr = element(patch,:);
    Fvalues(:,patch) = reshape(F,Nxi*Neta,1);
    Kvalues(:,patch) = reshape(Kloc,noDofsPatch^2,1);
    Findices(:,patch) = sctrGamma;
    spIdxRow(:,patch) = copyVector(sctr,noDofsPatch,1);
    spIdxCol(:,patch) = copyVector(sctr,noDofsPatch,2);
    spIdxM(:,patch) = sctr;
end
spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Kvalues = reshape(Kvalues,numel(Kvalues),1);
Mvalues = reshape(Mvalues,numel(Mvalues),1);
spIdx = [spIdxRow, spIdxCol];
clear spIdxRow spIdxCol
[spIdx,~,I] = unique(spIdx,'rows','stable');
Kvalues = accumarray(I,Kvalues);
noDofs_new = noDofs + noSurfDofs*(N-1);
varCol.A_K = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Kvalues,noDofs_new,noDofs_new,numel(Kvalues));
clear Kvalues

spIdxM = reshape(spIdxM,numel(spIdxM),1);
spIdx = [spIdxM, spIdxM];
[spIdx,~,I] = unique(spIdx,'rows','stable');
Mvalues = accumarray(I,Mvalues);

varCol.A_M = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Mvalues,noDofs_new,noDofs_new,numel(Mvalues));

varCol.FF = zeros(noDofs_new,no_angles);        % external force vector
for alpha_s_Nr = 1:no_angles
    varCol.FF(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),Findices,noDofs_new);
end
clear Fvalues

spIdxRowInf = reshape(spIdxRowInf,numel(spIdxRowInf),1);
spIdxColInf = reshape(spIdxColInf,numel(spIdxColInf),1);
Kinfvalues = reshape(Kinfvalues,numel(Kinfvalues),1);
spIdx = [spIdxRowInf, spIdxColInf];
newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRowInf));
clear spIdxRowInf spIdxColInf
[spIdx,~,I] = unique(spIdx,'rows','stable');
Kinfvalues = accumarray(I,Kinfvalues);
varCol.Ainf = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Kinfvalues,noDofs_new,noDofs_new,numel(Kinfvalues));
% [A_gamma_a, newDofsToRemove] = addInfElements3_testSEM(varCol);

varCol.dofsToRemove = sort(unique([varCol.dofsToRemove newDofsToRemove]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = delta(i,j)

if i == j
    d = 1;
else
    d = 0;
end

function invA = inv2(A,J)
% Formula on http://mathworld.wolfram.com/MatrixInverse.html
invA = zeros(size(A));
invA(1,1,:,:,:) =  A(2,2,:,:,:)./J;
invA(2,1,:,:,:) = -A(2,1,:,:,:)./J;
invA(1,2,:,:,:) = -A(1,2,:,:,:)./J;
invA(2,2,:,:,:) =  A(1,1,:,:,:)./J;


function A = multiplyByWeights(A,rho,Nxi,Neta,Nzeta,d)

A = A.*repmat(reshape(rho{1},1,1,Nxi,1,1),  d,d,1,  Neta,Nzeta);
A = A.*repmat(reshape(rho{2},1,1,1,Neta,1), d,d,Nxi,1,   Nzeta);
A = A.*repmat(reshape(rho{3},1,1,1,1,Nzeta),d,d,Nxi,Neta,1);

function A = multiplyByWeightsAtBoundary(A,rho,Nxi,Neta)

A = A.*repmat(reshape(rho{1},Nxi,1),  1,  Neta);
A = A.*repmat(reshape(rho{2},1,Neta), Nxi,1);

function A = multiplyByWeightsAtBoundary2(A,rho,Nxi,Neta,d1,d2)

A = A.*repmat(reshape(rho{1},1,1,Nxi,1),  d1, d2, 1,  Neta);
A = A.*repmat(reshape(rho{2},1,1,1,Neta), d1, d2, Nxi,1);