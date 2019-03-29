function [A, A2, A3, newDofsToRemove] = addInfElements3_ROM(varCol)

elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
patches = varCol.patches;
knotVecs = varCol.knotVecs;
noPatches = varCol.noPatches;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);

gluedNodes = varCol.gluedNodes;
N = varCol.N;
formulation = varCol.formulation;
A_2 = varCol.A_2;
x_0 = varCol.x_0;
noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

k = varCol.k(1);
Upsilon = varCol.Upsilon;
r_a = varCol.r_a;

D = varCol.D;
Dt = varCol.Dt;

noParams = 2;
noSurfDofs = 0;
noElemsPatch = zeros(noPatches,1);
noEl = zeros(noPatches,noParams);
for i = 1:noPatches
    n_xi = patches{i}.nurbs.number(1);
    n_eta = patches{i}.nurbs.number(2);
    noSurfDofs = noSurfDofs + n_xi*n_eta;
    for j = 1:noParams
        noEl(i,j) = size(patches{i}.elRange{j},1);
    end
    noElemsPatch(i) = size(patches{i}.elRange{1},1)*size(patches{i}.elRange{2},1);
end

zeta1Nodes = zeros(1,noSurfDofs);
counter = 1;
shiftIdx = 0;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    for j = 1:n_eta
        for i = 1:n_xi
            zeta1Nodes(counter) = shiftIdx + (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
            counter = counter + 1;
        end
    end
    shiftIdx = shiftIdx + n_zeta*n_eta*n_xi;
end
noElems = sum(noElemsPatch);
pIndex = zeros(noElems,1);
element = zeros(noElems,n_en);
index = zeros(noElems,noParams);
e = 1;
maxDof = 0;
jEl = zeros(1,2);
for i = 1:noPatches
    Xi = knotVecs{i}{1};
    Eta = knotVecs{i}{2};
    n_xi = patches{i}.nurbs.number(1);
    n_eta = patches{i}.nurbs.number(2);
    [surfElement, indexXiEta, noElemsXiEta] = generateIGA2DMesh(Xi, Eta, p_xi, p_eta, n_xi, n_eta);
    index(e:e+noElemsPatch(i)-1,:) = indexXiEta + repmat(jEl,noElemsXiEta,1);    
    pIndex(e:e+noElemsXiEta-1) = i;
    element(e:e+noElemsXiEta-1,:) = maxDof + surfElement;
    jEl = jEl + noEl(i,:);
    maxDof = maxDof + n_xi*n_eta;
    e = e + noElemsPatch(i);
end
element2 = element;
% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    indices = (zeta1Nodes(element(:)) == gluedNodes{i}(1));
    parentIdx = element(indices);
    element(indices) = parentIdx;
    for j = 2:length(gluedNodes{i})
        indices = (zeta1Nodes(element(:)) == gluedNodes{i}(j));
        element(indices) = parentIdx;
    end
end


%% Evaluate analytic integrals in ``radial'' direction. 
% Note that the last two integrals (I1(end) and I1(end-1),
% I2(end) and I2(end-1)) will be redundant for the cases 'BGC' and 'BGU'
B1 = zeros(2*N+4,1);
B2 = zeros(2*N+3,1);
varrho1 = Upsilon/r_a;
varrho2 = k*r_a;
varrho3 = k*Upsilon;
for n = 1:2*N+4
    B1(n) = radialIntegral3(n, varrho1, varrho2, formulation, 1);
    if n < 2*N+4
        B2(n) = radialIntegral3(n, varrho1, varrho2, formulation, 2);
    end
end

%% Calculate contribution from infinite elements

spIdxRow = zeros((N*n_en)^2,noElems);
spIdxCol = zeros((N*n_en)^2,noElems);
Avalues = zeros((N*n_en)^2,noElems);
A2values = zeros((N*n_en)^2,noElems);
A3values = zeros((N*n_en)^2,noElems);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
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

    sctrLocal = element(e,:);
    sctrGlobal = zeta1Nodes(sctrLocal);

    pts = controlPts(sctrGlobal,:);
    wgts = weights(zeta1Nodes(element2(e,:)),:); % New
    
    A1_IJ = zeros(n_en);
    A2_IJ = zeros(n_en);
    A4_IJ = zeros(n_en);
    A5_IJ = zeros(n_en);
    A3_IJ = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        J = pts'*[dRdxi' dRdeta'];
        
        X = R*pts;
        Xt = A_2*(X-x_0)';
        xt = Xt(1);
        yt = Xt(2);
        zt = Xt(3);
        
        [~, theta, ~, d1, d2] = evaluateProlateCoords(xt,yt,zt,Upsilon);
        
        if sin(theta) < 10*eps
            temp = R'*R*(J(1,1)*J(2,2)-J(1,2)*J(2,1))/(r_a^2-Upsilon^2);
            tempI = ( (J(1,1)^2+J(2,1)^2)*(dRdeta'*dRdeta) ...
                     -(J(1,1)*J(1,2)+J(2,1)*J(2,2))*(dRdeta'*dRdxi + dRdxi'*dRdeta) ...
                     +(J(1,2)^2+J(2,2)^2)*(dRdxi'*dRdxi))/(J(1,1)*J(2,2)-J(1,2)*J(2,1));
        end
        if theta < 10*eps
            A1_IJ = A1_IJ + temp*J_2*wt; 
            A2_IJ = A2_IJ + tempI*J_2*wt;   
            A3_IJ = A3_IJ + temp*J_2*wt; % note that cos(theta)^2 = 1 in this case
        elseif theta > pi-10*eps
            A1_IJ = A1_IJ - temp*J_2*wt; 
            A2_IJ = A2_IJ - tempI*J_2*wt;   
            A3_IJ = A3_IJ - temp*J_2*wt; % note that cos(theta)^2 = 1 in this case
        else
            DPDX = dPdX(xt,yt,zt,Upsilon,r_a,d1,d2)*A_2;

            J3 = DPDX(2:3,:)*J(:,1:2);
            J_3 = abs(det(J3));

            dRdP = J3'\[dRdxi; dRdeta];
            dRdtheta = dRdP(1,:);
            dRdphi = dRdP(2,:);
            
            A1_IJ = A1_IJ + R'*R*sin(theta)                               *J_3*J_2*wt;  
            A2_IJ = A2_IJ + dRdtheta'*dRdtheta*sin(theta)                 *J_3*J_2*wt;  
            A3_IJ = A3_IJ + R'*R*cos(theta)^2*sin(theta)                  *J_3*J_2*wt;  
            A4_IJ = A4_IJ + dRdphi'*dRdphi/sin(theta)                     *J_3*J_2*wt;  
            A5_IJ = A5_IJ + dRdphi'*dRdphi*cos(theta)^2/sin(theta)        *J_3*J_2*wt;  
        end
    end   
    
    A_inf_values_temp = zeros((N*n_en)^2,1);
    A2_inf_values_temp = zeros((N*n_en)^2,1);
    A3_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    for m = 1:N
        for n = 1:N
            temp = zeros(n_en);
            temp2 = zeros(n_en);
            temp3 = zeros(n_en);
            for nt = 1:N
                for mt = 1:N
                    switch formulation
                        case 'PGU'
                            tempI  = A1_IJ*(-2*varrho2^2*B1(nt+mt) - 1i*varrho2*(nt+mt+2)*B1(nt+mt+1) + ((nt+2)*mt + varrho3^2)*B1(nt+mt+2) ...
                                          +1i*varrho1*varrho3*(nt+mt+2)*B1(nt+mt+3) - varrho1^2*(nt+2)*mt*B1(nt+mt+4)) ...
                                   + A2_IJ*B1(nt+mt+2) + varrho3^2*A3_IJ*B1(nt+mt+2) ...
                                   + A4_IJ*B2(nt+mt+1) - varrho1^2*A5_IJ*B2(nt+mt+3);
                        case 'PGC'
                            tempI  = A1_IJ*(- 1i*varrho2*(nt-mt+2)*B1(nt+mt+1) + ((nt+2)*mt - varrho3^2)*B1(nt+mt+2) ...
                                          +1i*varrho1*varrho3*(nt-mt+2)*B1(nt+mt+3) - varrho1^2*(nt+2)*mt*B1(nt+mt+4)) ...
                                   + A2_IJ*B1(nt+mt+2) + varrho3^2*A3_IJ*B1(nt+mt+2) ...
                                   + A4_IJ*B2(nt+mt+1) - varrho1^2*A5_IJ*B2(nt+mt+3);
                        case 'BGU'
                            if mt+nt == 2
                                tempI  = A1_IJ*(- 2*1i*varrho2*B1(1) + (1 + varrho3^2)*B1(2) ...
                                              +2*1i*varrho1*varrho3*B1(3) - varrho1^2*B1(4) - 1i*varrho2*exp(2*1i*varrho2)) ...
                                       + A2_IJ*B1(2) + varrho3^2*A3_IJ*B1(2) ...
                                       + A4_IJ*B2(1) - varrho1^2*A5_IJ*B2(3);
                            else
                                tempI  = A1_IJ*(-2*varrho2^2*B1(nt+mt-2) - 1i*varrho2*(nt+mt)*B1(nt+mt-1) + (nt*mt + varrho3^2)*B1(nt+mt) ...
                                              +1i*varrho1*varrho3*(nt+mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2)) ...
                                       + A2_IJ*B1(nt+mt) + varrho3^2*A3_IJ*B1(nt+mt) ...
                                       + A4_IJ*B2(nt+mt-1) - varrho1^2*A5_IJ*B2(nt+mt+1);
                            end
                        case 'BGC'
                            if mt+nt == 2
                                tempI  = A1_IJ*((1 - varrho3^2)*B1(2) ...
                                              - varrho1^2*B1(4) - 1i*varrho2) ...
                                       + A2_IJ*B1(2) + varrho3^2*A3_IJ*B1(2) ...
                                       + A4_IJ*B2(1) - varrho1^2*A5_IJ*B2(3);
                                tempI2  = A1_IJ*( - 2*Upsilon^2*k*B1(2) - 1i*r_a) ...
                                       + 2*Upsilon^2*k*A3_IJ*B1(2);
                                tempI3  = 2*Upsilon^2*(A3_IJ-A1_IJ)*B1(2);
                            else
                                tempI  = A1_IJ*(-1i*varrho2*(nt-mt)*B1(nt+mt-1) + (nt*mt - varrho3^2)*B1(nt+mt) ...
                                              +1i*varrho1*varrho3*(nt-mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2)) ...
                                       + A2_IJ*B1(nt+mt) + varrho3^2*A3_IJ*B1(nt+mt) ...
                                       + A4_IJ*B2(nt+mt-1) - varrho1^2*A5_IJ*B2(nt+mt+1);
                                tempI2  = A1_IJ*(-1i*r_a*(nt-mt)*B1(nt+mt-1) - 2*Upsilon^2*k*B1(nt+mt) ...
                                              +1i*varrho1*Upsilon*(nt-mt)*B1(nt+mt+1)) ...
                                       + 2*Upsilon^2*k*A3_IJ*B1(nt+mt);
                                tempI3  = 2*Upsilon^2*(A3_IJ-A1_IJ)*B1(nt+mt);
                            end
                    end
                    
                    switch formulation
                        case {'PGU', 'BGU'}
                            temp = temp + Dt(n,nt)*D(m,mt)*r_a*tempI*exp(-2*1i*varrho2);
                            temp2 = temp2 + Dt(n,nt)*D(m,mt)*r_a*tempI*exp(-2*1i*varrho2);
                            temp3 = temp3 + Dt(n,nt)*D(m,mt)*r_a*tempI*exp(-2*1i*varrho2);
                        case {'PGC', 'BGC'}
                            temp = temp + Dt(n,nt)*D(m,mt)*r_a*tempI;
                            temp2 = temp2 + Dt(n,nt)*D(m,mt)*r_a*tempI2;
                            temp3 = temp3 + Dt(n,nt)*D(m,mt)*r_a*tempI3;
                    end
                end
            end
            indices = (1:n_en^2)+(n_en^2*N*(m-1) + n_en^2*(n-1));
            A_inf_values_temp(indices) = reshape(temp,n_en^2,1);
            A2_inf_values_temp(indices) = reshape(temp2,n_en^2,1);
            A3_inf_values_temp(indices) = reshape(temp3,n_en^2,1);
            if n == 1
                IEsctr = sctrGlobal;
            else
                IEsctr = noDofs+sctrLocal+noSurfDofs*(n-2);
            end
            spIdxRow_temp(indices) = copyVector(IEsctr,n_en,1);
            if m == 1
                IEsctr = sctrGlobal;
            else
                IEsctr = noDofs+sctrLocal+noSurfDofs*(m-2);
            end
            spIdxCol_temp(indices) = copyVector(IEsctr,n_en,2);
        end
    end
    Avalues(:,e) = A_inf_values_temp;
    A2values(:,e) = A2_inf_values_temp;
    A3values(:,e) = A3_inf_values_temp;
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noSurfDofs*(N-1);

spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Avalues = reshape(Avalues,numel(Avalues),1);
A2values = reshape(A2values,numel(A2values),1);
A3values = reshape(A3values,numel(A3values),1);

[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
Avalues = accumarray(IuniqueIdx,Avalues);
A2values = accumarray(IuniqueIdx,A2values);
A3values = accumarray(IuniqueIdx,A3values);

A = sparse(spIdx(:,1),spIdx(:,2),Avalues,noDofs_new,noDofs_new,numel(IuniqueIdx));
A2 = sparse(spIdx(:,1),spIdx(:,2),A2values,noDofs_new,noDofs_new,numel(IuniqueIdx));
A3 = sparse(spIdx(:,1),spIdx(:,2),A3values,noDofs_new,noDofs_new,numel(IuniqueIdx));

newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));

