function [A, newDofsToRemove] = addInfElements3(varCol)

elRange = varCol.elRange(1:2);
knotVecs = varCol.knotVecs;

degree = varCol.degree(1:2); % assume degree is equal in all patches

N = varCol.N;
formulation = varCol.formulation;
A_2 = varCol.A_2;
x_0 = varCol.x_0;
noDofs = varCol.noDofs;

d_p = varCol.patches{1}.nurbs.d_p;

weights = varCol.weights;
controlPts = varCol.controlPts;

k = varCol.k(1);
Upsilon = varCol.Upsilon;
r_a = varCol.r_a;

D = varCol.D;
Dt = varCol.Dt;
[nodes, noElems, element, element2, index, pIndex, n_en, noSurfDofs] = meshBoundary(varCol,1);


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

[Q, W] = gaussTensorQuad(degree+1);
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(d_p-1,2);
    for i = 1:d_p-1
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^(d_p-1);
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    
    
    sctrLocal = element(e,:);
    sctrGlobal = nodes(sctrLocal);

    pts = controlPts(sctrGlobal,:);
    wgts = weights(nodes(element2(e,:)),:); % New
    
    A1_IJ = zeros(n_en);
    A2_IJ = zeros(n_en);
    A4_IJ = zeros(n_en);
    A5_IJ = zeros(n_en);
    A3_IJ = zeros(n_en);
            
    R = NURBSbasis(I, xi, degree, knots, wgts);
    dXdxi = R{2}*pts;
    dXdeta = R{3}*pts;
    a11 = dXdxi(:,1);
    a21 = dXdxi(:,2);
    a31 = dXdxi(:,3);
    a12 = dXdeta(:,1);
    a22 = dXdeta(:,2);
    a32 = dXdeta(:,3);

    X = R{1}*pts;
    Xt = (X-x_0)*A_2.';

    [~, theta_arr, ~, d1, d2] = evaluateProlateCoords(Xt,Upsilon);
        
    for gp = 1:numel(W)
        theta = theta_arr(gp);
        RR = R{1}(gp,:)'*R{1}(gp,:);
        if sin(theta) < 10*eps
            temp = RR*(a11(gp)*a22(gp)-a12(gp)*a21(gp))/(r_a^2-Upsilon^2);
            temp2 = ( (a11(gp)^2+a21(gp)^2)*(R{3}(gp,:)'*R{3}(gp,:)) ...
                     -(a11(gp)*a12(gp)+a21(gp)*a22(gp))*(R{3}(gp,:)'*R{2}(gp,:) + R{2}(gp,:)'*R{3}(gp,:)) ...
                     +(a12(gp)^2+a22(gp)^2)*(R{2}(gp,:)'*R{2}(gp,:)))/(a11(gp)*a22(gp)-a12(gp)*a21(gp));
        end
        if theta < 10*eps
            A1_IJ = A1_IJ + temp*J_2*W(gp); 
            A2_IJ = A2_IJ + temp2*J_2*W(gp);   
            A3_IJ = A3_IJ + temp*J_2*W(gp); % note that cos(theta)^2 = 1 in this case
        elseif theta > pi-10*eps
            A1_IJ = A1_IJ - temp*J_2*W(gp); 
            A2_IJ = A2_IJ - temp2*J_2*W(gp);   
            A3_IJ = A3_IJ - temp*J_2*W(gp); % note that cos(theta)^2 = 1 in this case
        else
            DPDX = dPdX(Xt(gp,:),Upsilon,r_a,d1(gp),d2(gp))*A_2;
            J = [a11(gp) a12(gp);
                 a21(gp) a22(gp);
                 a31(gp) a32(gp)];
                 
            J3 = DPDX(2:3,:)*J;
            J_3 = abs(det(J3));

            dRdP = J3'\[R{2}(gp,:); R{3}(gp,:)];
            dRdtheta = dRdP(1,:);
            dRdphi = dRdP(2,:);
            
            A1_IJ = A1_IJ + RR*sin(theta)                          	*J_3*J_2*W(gp);  
            A2_IJ = A2_IJ + dRdtheta'*dRdtheta*sin(theta)        	*J_3*J_2*W(gp);  
            A3_IJ = A3_IJ + RR*cos(theta)^2*sin(theta)           	*J_3*J_2*W(gp);  
            A4_IJ = A4_IJ + dRdphi'*dRdphi/sin(theta)              	*J_3*J_2*W(gp);  
            A5_IJ = A5_IJ + dRdphi'*dRdphi*cos(theta)^2/sin(theta) 	*J_3*J_2*W(gp);  
        end
    end   
    
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    for m = 1:N
        for n = 1:N
            temp = zeros(n_en);
            for nt = 1:N
                for mt = 1:N
                    switch formulation
                        case 'PGU'
                            temp2  = A1_IJ*(-2*varrho2^2*B1(nt+mt) - 1i*varrho2*(nt+mt+2)*B1(nt+mt+1) + ((nt+2)*mt + varrho3^2)*B1(nt+mt+2) ...
                                          +1i*varrho1*varrho3*(nt+mt+2)*B1(nt+mt+3) - varrho1^2*(nt+2)*mt*B1(nt+mt+4)) ...
                                   + A2_IJ*B1(nt+mt+2) + varrho3^2*A3_IJ*B1(nt+mt+2) ...
                                   + A4_IJ*B2(nt+mt+1) - varrho1^2*A5_IJ*B2(nt+mt+3);
                        case 'PGC'
                            temp2  = A1_IJ*(- 1i*varrho2*(nt-mt+2)*B1(nt+mt+1) + ((nt+2)*mt - varrho3^2)*B1(nt+mt+2) ...
                                          +1i*varrho1*varrho3*(nt-mt+2)*B1(nt+mt+3) - varrho1^2*(nt+2)*mt*B1(nt+mt+4)) ...
                                   + A2_IJ*B1(nt+mt+2) + varrho3^2*A3_IJ*B1(nt+mt+2) ...
                                   + A4_IJ*B2(nt+mt+1) - varrho1^2*A5_IJ*B2(nt+mt+3);
                        case 'BGU'
                            if mt+nt == 2
                                temp2  = A1_IJ*(- 2*1i*varrho2*B1(1) + (1 + varrho3^2)*B1(2) ...
                                              +2*1i*varrho1*varrho3*B1(3) - varrho1^2*B1(4) - 1i*varrho2*exp(2*1i*varrho2)) ...
                                       + A2_IJ*B1(2) + varrho3^2*A3_IJ*B1(2) ...
                                       + A4_IJ*B2(1) - varrho1^2*A5_IJ*B2(3);
                            else
                                temp2  = A1_IJ*(-2*varrho2^2*B1(nt+mt-2) - 1i*varrho2*(nt+mt)*B1(nt+mt-1) + (nt*mt + varrho3^2)*B1(nt+mt) ...
                                              +1i*varrho1*varrho3*(nt+mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2)) ...
                                       + A2_IJ*B1(nt+mt) + varrho3^2*A3_IJ*B1(nt+mt) ...
                                       + A4_IJ*B2(nt+mt-1) - varrho1^2*A5_IJ*B2(nt+mt+1);
                            end
                        case 'BGC'
                            if mt+nt == 2
                                temp2  = A1_IJ*((1 - varrho3^2)*B1(2) ...
                                              - varrho1^2*B1(4) - 1i*varrho2) ...
                                       + A2_IJ*B1(2) + varrho3^2*A3_IJ*B1(2) ...
                                       + A4_IJ*B2(1) - varrho1^2*A5_IJ*B2(3);
                            else
                                temp2  = A1_IJ*(-1i*varrho2*(nt-mt)*B1(nt+mt-1) + (nt*mt - varrho3^2)*B1(nt+mt) ...
                                              +1i*varrho1*varrho3*(nt-mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2)) ...
                                       + A2_IJ*B1(nt+mt) + varrho3^2*A3_IJ*B1(nt+mt) ...
                                       + A4_IJ*B2(nt+mt-1) - varrho1^2*A5_IJ*B2(nt+mt+1);
                            end
                    end
                    
                    switch formulation
                        case {'PGU', 'BGU'}
                            temp = temp + Dt(n,nt)*D(m,mt)*r_a*temp2*exp(-2*1i*varrho2);
                        case {'PGC', 'BGC'}
                            temp = temp + Dt(n,nt)*D(m,mt)*r_a*temp2;
                    end
                end
            end
            indices = (1:n_en^2)+(n_en^2*N*(m-1) + n_en^2*(n-1));
            A_inf_values_temp(indices) = reshape(temp,n_en^2,1);
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
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noSurfDofs*(N-1);

spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
Avalues = accumarray(IuniqueIdx,Avalues);

A = sparse(spIdx(:,1),spIdx(:,2),Avalues,noDofs_new,noDofs_new,numel(IuniqueIdx));

newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));

