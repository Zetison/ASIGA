function [A, newDofsToRemove] = infElementsNonSepGeom(varCol)

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

extraGP = varCol.extraGP;

N = varCol.N;
formulation = varCol.formulation;

noDofs = varCol.noDofs;


%Use chebychev polynomials
chimin = varCol.chimin;
chimax = varCol.chimax;
x_0 = varCol.x_0;
A_2 = varCol.A_2;
D = varCol.D;
Dt = varCol.Dt;

k = varCol.k;
Upsilon = varCol.Upsilon;

useApproxRadialIntegrals = 2;
B1splines = {};
B2splines = {};
B1cheb = {};
B2cheb = {};

if useApproxRadialIntegrals == 1
    p_splineApprox = 2;
    B1splines = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, formulation, 1),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
    B2splines = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, formulation, 2),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
elseif useApproxRadialIntegrals == 2
    B1cheb = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                formulation, 1), 1e-6);
    end
    B2cheb = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                formulation, 2), 1e-6);
    end
end

%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

spIdxRow = zeros((N*n_en)^2,noElems);
spIdxCol = zeros((N*n_en)^2,noElems);
Avalues = zeros((N*n_en)^2,noElems);
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);

% max_r_a_recorded = -inf;
% min_r_a_recorded = inf;

% totArea = 0;
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:)); % New   
    idXi = index(e,1); 
    idEta = index(e,2);

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    temp = zeros(n_en,n_en,N,N);
    for gp = 1:size(W2D,1)
%         tic
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        
        J = pts.'*[dRdxi.' dRdeta.'];
        
        X = R*pts;
        Xt = A_2*(X-x_0).';
        xt = Xt(1);
        yt = Xt(2);
        zt = Xt(3);
        
        [r_a, theta, ~, c1, c2] = evaluateProlateCoords(xt,yt,zt,Upsilon);
        
        if r_a > chimax || r_a < chimin
            error('chi is not in range [chimin chimax]')
        end
%         if r_a < min_r_a_recorded
%             min_r_a_recorded = r_a;
%         end
%         if r_a > max_r_a_recorded
%             max_r_a_recorded = r_a;
%         end
        DPDX = dPdX(xt,yt,zt,Upsilon,r_a,c1,c2)*A_2;
        
        J3 = DPDX(2:3,:)*J;
        J_3 = abs(det(J3));
        
        
        [dr_adtheta, dr_adphi] = dchidP(DPDX(1,:),J,J3);
        
        dRdP = J3'\[dRdxi; dRdeta];
        
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);
        if useApproxRadialIntegrals == 1
            B1 = zeros(2*N+4,1);
            xi_pt = (r_a-chimin)/(chimax-chimin);
            for i = 1:2*N+4
                B1(i) = evaluateNURBS(B1splines{i}, xi_pt);
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                B2(i) = evaluateNURBS(B2splines{i}, xi_pt);
            end
        elseif useApproxRadialIntegrals == 2
            B1 = zeros(2*N+4,1);
            x_cheb = (2*r_a-chimax-chimin)/(chimax-chimin);
            
            for i = 1:2*N+4
                for j = 1:length(B1cheb{i})
                    if j == 1 || j == length(B1cheb{i})
                        B1(i) = B1(i) + 0.5*B1cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                    else
                        B1(i) = B1(i) + B1cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                    end
                end
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                for j = 1:length(B2cheb{i})
                    if j == 1 || j == length(B2cheb{i})
                        B2(i) = B2(i) + 0.5*B2cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                    else
                        B2(i) = B2(i) + B2cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                    end
                end
            end      
        else
            B1 = zeros(2*N+4,1);
            for i = 1:2*N+4
                B1(i) = radialIntegral2(i, r_a, k, Upsilon, formulation, 1);
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                B2(i) = radialIntegral2(i, r_a, k, Upsilon, formulation, 2);
            end      
        end

        RR = R'*R;
        RtRt = dRdtheta'*dRdtheta;  
        RpRp = dRdphi'*dRdphi;  
        
        RtR = dRdtheta'*R;  
        RRt = R'*dRdtheta; 
        RpR = dRdphi'*R; 
        RRp = R'*dRdphi;

        varrho1 = Upsilon/r_a;
        varrho2 = k*r_a;
        varrho3 = k*Upsilon;
        for nt = 1:N
            for mt = 1:N                
                switch formulation
                    case 'PGU'
                        temp2 = RR*(-2*varrho2^2*B1(nt+mt) - 1i*varrho2*(nt+mt+2)*B1(nt+mt+1) +((nt+2)*mt + varrho3^2*(1+cos(theta)^2))*B1(nt+mt+2) ...
                                    +1i*varrho1*varrho3*(nt+mt+2)*B1(nt+mt+3)-varrho1^2*(nt+2)*mt*B1(nt+mt+4) ...
                                    +(-varrho2^2-1i*varrho2*(nt+mt+2)+(nt+2)*mt) ...
                                        *1/r_a^2*(dr_adtheta^2*B1(nt+mt+2)+dr_adphi^2/sin(theta)^2*(B2(nt+mt+1)-varrho1^2*cos(theta)^2*B2(nt+mt+3)))) ...
                                +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt+2-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(nt+mt+2) ...
                                +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt+2-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                    *1/sin(theta)^2*(B2(nt+mt+1) - varrho1^2*cos(theta)^2*B2(nt+mt+3));
                    case 'PGC'
                        temp2 = RR*(- 1i*varrho2*(nt-mt+2)*B1(nt+mt+1) + ((nt+2)*mt + varrho3^2*(-1+cos(theta)^2))*B1(nt+mt+2) ...
                                    +1i*varrho1*varrho3*(nt-mt+2)*B1(nt+mt+3)-varrho1^2*(nt+2)*mt*B1(nt+mt+4) ...
                                    +(varrho2^2-1i*varrho2*(nt-mt+2)+(nt+2)*mt) ...
                                        *1/r_a^2*(dr_adtheta^2*B1(nt+mt+2)+dr_adphi^2/sin(theta)^2*(B2(nt+mt+1)-varrho1^2*cos(theta)^2*B2(nt+mt+3)))) ...
                                +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt+2+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(nt+mt+2) ...
                                +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt+2+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                    *1/sin(theta)^2*(B2(nt+mt+1) - varrho1^2*cos(theta)^2*B2(nt+mt+3));
                    case 'BGU'
                        if mt+nt == 2
                            temp2 = RR*(-2*1i*varrho2*B1(1) +(1 + varrho3^2*(1+cos(theta)^2))*B1(2) ...
                                        +2*1i*varrho1*varrho3*B1(3)-varrho1^2*B1(4) - 1i*varrho2*exp(2*1i*varrho2) ...
                                        +(-varrho2^2-2*1i*varrho2+1) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(2)+dr_adphi^2/sin(theta)^2*(B2(1)-varrho1^2*cos(theta)^2*B2(3)))) ...
                                    +(RtRt + (1-1i*varrho2)/r_a*dr_adtheta*RtR + (1-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(2) ...
                                    +(RpRp + (1-1i*varrho2)/r_a*dr_adphi*RpR + (1-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(1) - varrho1^2*cos(theta)^2*B2(3));
                        else
                            temp2 = RR*(-2*varrho2^2*B1(nt+mt-2) - 1i*varrho2*(nt+mt)*B1(nt+mt-1) +(nt*mt + varrho3^2*(1+cos(theta)^2))*B1(nt+mt) ...
                                        +1i*varrho1*varrho3*(nt+mt)*B1(nt+mt+1)-varrho1^2*nt*mt*B1(nt+mt+2) ...
                                        +(-varrho2^2-1i*varrho2*(nt+mt)+nt*mt) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(nt+mt)+dr_adphi^2/sin(theta)^2*(B2(nt+mt-1)-varrho1^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                    +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(nt+mt) ...
                                    +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(nt+mt-1) - varrho1^2*cos(theta)^2*B2(nt+mt+1));
                        end
                    case 'BGC'
                        if mt+nt == 2
                            temp2 = RR*((1+varrho3^2*(-1+cos(theta)^2))*B1(2) ...
                                        -varrho1^2*B1(4) - 1i*varrho2...
                                        +(varrho2^2+1) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(2)+dr_adphi^2/sin(theta)^2*(B2(1)-varrho1^2*cos(theta)^2*B2(3)))) ...
                                    +(RtRt + (1-1i*varrho2)/r_a*dr_adtheta*RtR+(1+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(2) ...
                                    +(RpRp + (1-1i*varrho2)/r_a*dr_adphi*RpR+(1+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(1) - varrho1^2*cos(theta)^2*B2(3));
                        else
                            temp2 = RR*(-1i*varrho2*(nt-mt)*B1(nt+mt-1) + (nt*mt+varrho3^2*(-1+cos(theta)^2))*B1(nt+mt) ...
                                        +1i*varrho1*varrho3*(nt-mt)*B1(nt+mt+1)-varrho1^2*nt*mt*B1(nt+mt+2) ...
                                        +(varrho2^2-1i*varrho2*(nt-mt)+nt*mt) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(nt+mt)+dr_adphi^2/sin(theta)^2*(B2(nt+mt-1)-varrho1^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                    +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR+(nt+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(nt+mt) ...
                                    +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR+(nt+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(nt+mt-1) - varrho1^2*cos(theta)^2*B2(nt+mt+1));
                        end
                end
                for n = 1:N
                    for m = 1:N
                        if Dt(n,nt)*D(m,mt) ~= 0
                            switch formulation
                                case {'PGU', 'BGU'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*exp(-2*1i*varrho2)*sin(theta)*J_3*J_2*wt;
                                case {'PGC', 'BGC'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*sin(theta)*J_3*J_2*wt;
                            end
                        end
                    end
                end
            end
        end
%         totArea = totArea + J_3*J_2*wt;
    end  
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    counter = 1;
    for n = 1:N
        for m = 1:N
            indices = counter:(counter+n_en^2-1);
            A_inf_values_temp(indices) = reshape(temp(:,:,n,m),n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctr+(noDofs*(n-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctr+(noDofs*(m-1)),n_en,2);
            counter = counter + n_en^2;
        end
    end
    
    Avalues(:,e) = A_inf_values_temp;
    
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end
% errorInTotArea = totArea-2*pi^2;

% max_r_a_recorded
% min_r_a_recorded
noDofs_new = noDofs*N;

spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
Avalues = accumarray(IuniqueIdx,Avalues);

A = sparse(spIdx(:,1),spIdx(:,2),Avalues,noDofs_new,noDofs_new,numel(IuniqueIdx));

if min(size(A)) < noDofs_new
    A(noDofs_new, noDofs_new) = 0;
end
newDofsToRemove = setdiff(1:noDofs_new,unique(spIdxRow));

