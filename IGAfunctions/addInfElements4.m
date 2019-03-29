function [A_inf, newDofsToRemove] = addInfElements4(varCol, k, Upsilon)

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

n_xi = varCol.nurbs.number(1);
n_eta = varCol.nurbs.number(2);
n_zeta = varCol.nurbs.number(3);
gluedNodes = varCol.gluedNodes;
N = varCol.N;
formulation = varCol.formulation;
% rm = varCol.rm;

noSurfDofs = n_xi*n_eta;
noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;

zeta1Nodes = zeros(1,n_xi*n_eta);
counter = 1;
for j = 1:n_eta
    for i = 1:n_xi
        zeta1Nodes(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
        counter = counter + 1;
    end
end

[XiEtaMesh, indexXiEta, noElemsXiEta] = generateIGA2DMesh(Xi, Eta, p_xi, p_eta, n_xi, n_eta);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (XiEtaMesh == gluedNodes{i}(j));
        XiEtaMesh(indices) = parentIdx;
    end
end

%% Find coefficients for radial shape functions in infinite elements
ramin = varCol.chimin;
ramax = varCol.chimax;
x_0 = varCol.x_0;
A_2 = varCol.A_2;
D = varCol.D;
Dt = varCol.Dt;


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
                                           ramin,ramax,max(ceil((ramax-ramin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
    B2splines = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, formulation, 2),...
                                           ramin,ramax,max(ceil((ramax-ramin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
elseif useApproxRadialIntegrals == 2
    B1cheb = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, ramin+0.5*(ramax-ramin)*(1+x), k, Upsilon, ...
                                                formulation, 1), 1e-6);
    end
    B2cheb = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, ramin+0.5*(ramax-ramin)*(1+x), k, Upsilon, ...
                                                formulation, 2), 1e-6);
    end
end

%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

spIdxRow = zeros((N*n_en)^2,noElemsXiEta);
spIdxCol = zeros((N*n_en)^2,noElemsXiEta);
A_inf_values = zeros((N*n_en)^2,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
% for e = 1:noElemsXiEta
parfor e = 1:noElemsXiEta
    idXi = indexXiEta(e,1); 
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    
    temp = zeros(n_en,n_en,N,N);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights(zeta1Nodes));
        J = pts'*[dRdxi' dRdeta'];
        
        X = R*pts;
        Xt = A_2*(X-x_0)';
        xt = Xt(1);
        yt = Xt(2);
        zt = Xt(3);
        
        [r_a, theta, ~, c1, c2] = evaluateProlateCoords(xt,yt,zt,Upsilon);
        
%         if abs(chi-r) > 1e-13
%             keyboard
%         end
% chi
        if r_a > ramax || r_a < ramin
            error('r_a is not in range [ramin ramax]')
        end
%         if chi < minchi_recorded
%             minchi_recorded = chi;
%         end
%         if chi > maxchi_recorded,
%             maxchi_recorded = chi;
%         end
%         DPDX = dPdX(xt,yt,zt,Upsilon,r,theta,T);
        DPDX = dPdX(xt,yt,zt,Upsilon,r_a,c1,c2)*A_2;
        
        J3 = DPDX(2:3,:)*J;
        J_3 = abs(det(J3));
        
        
        [dradtheta, dradphi] = dchidP(DPDX(1,:),J,J3);
        
        dRdP = J3'\[dRdxi; dRdeta];
        
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);
%         toc
%         tic

        if useApproxRadialIntegrals == 1
            B1 = zeros(2*N+4,1);
            xi_pt = (r_a-ramin)/(ramax-ramin);
            for i = 1:2*N+4
                B1(i) = evaluateNURBS(B1splines{i}, xi_pt);
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                B2(i) = evaluateNURBS(B2splines{i}, xi_pt);
            end
        elseif useApproxRadialIntegrals == 2
            B1 = zeros(2*N+4,1);
            x_cheb = (2*r_a-ramax-ramin)/(ramax-ramin);
            
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
        
        for nt = 1:N
            for mt = 1:N                
                switch formulation
                    case 'PGU'
                        temp2 = RR*(-2*k^2*r_a^2*B1(nt+mt)-1i*k*r_a*(nt+mt+2)*B1(nt+mt+1)+(mt*(nt+2)+k^2*Upsilon^2)*B1(nt+mt+2) ...
                                   +1i*k*Upsilon^2/r_a*(nt+mt+2)*B1(nt+mt+3)-mt*(nt+2)*Upsilon^2/r_a^2*B1(nt+mt+4)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt+2) ...
                                   +1/r_a^2*(-k^2*r_a^2-1i*k*r_a*(nt+mt+2)+mt*(nt+2)) ...
                                           *(dradtheta^2*B1(nt+mt+2)+dradphi^2/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3)))) ...
                               +(RtRt+RtR*dradtheta*(-1i*k+mt/r_a)+RRt*dradtheta*(-1i*k+(nt+2)/r_a))*B1(nt+mt+2) ...
                               +RpRp/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3)) ...
                               +dradphi/sin(theta)^2*(RpR*(-1i*k+mt/r_a)+RRp*(-1i*k+(nt+2)/r_a)) ...
                                                     *(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3));

                    case 'PGC'
                        temp2 = RR*(-1i*k*r_a*(nt-mt+2)*B1(nt+mt+1)+(mt*(nt+2)-k^2*Upsilon^2)*B1(nt+mt+2) ...
                                   +1i*k*Upsilon^2/r_a*(nt-mt+2)*B1(nt+mt+3)-mt*(nt+2)*Upsilon^2/r_a^2*B1(nt+mt+4)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt+2) ...
                                   +1/r_a^2*(k^2*r_a^2-1i*k*r_a*(nt-mt+2)+mt*(nt+2)) ...
                                           *(dradtheta^2*B1(nt+mt+2)+dradphi^2/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3)))) ...
                               +(RtRt+RtR*dradtheta*(-1i*k+mt/r_a)+RRt*dradtheta*(1i*k+(nt+2)/r_a))*B1(nt+mt+2) ...
                               +RpRp/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3)) ...
                               +dradphi/sin(theta)^2*(RpR*(-1i*k+mt/r_a)+RRp*(1i*k+(nt+2)/r_a)) ...
                                                     *(B2(nt+mt+1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+3));
                    case 'BGU'
                        if mt+nt == 2
                            temp2 = RR*(-2*1i*k*r_a*B1(1)+(1+k^2*Upsilon^2)*B1(2) ...
                                        +2*1i*k*Upsilon^2/r_a*B1(3)-Upsilon^2/r_a^2*B1(4)+k^2*Upsilon^2*cos(theta)^2*B1(2) ...
                                        +1/r_a^2*(-k^2*r_a^2-2*1i*k*r_a+1) ...
                                           *(dradtheta^2*B1(2)+dradphi^2/sin(theta)^2*(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3))) ...
                                        -1i*k*exp(2*1i*k*r_a)*r_a) ...
                                    +(RtRt + RtR*dradtheta*(-1i*k+1/r_a)+RRt*dradtheta*(-1i*k+1/r_a))*B1(2) ...
                                    +RpRp/sin(theta)^2*(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3)) ...
                                    +dradphi/sin(theta)^2*(RpR*(-1i*k+1/r_a)+RRp*(-1i*k+1/r_a)) ...
                                                           *(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3));
                        else
                            temp2 = RR*(-2*k^2*r_a^2*B1(nt+mt-2)-1i*k*r_a*(nt+mt)*B1(nt+mt-1)+(mt*nt+k^2*Upsilon^2)*B1(nt+mt) ...
                                        +1i*k*Upsilon^2/r_a*(nt+mt)*B1(nt+mt+1)-mt*nt*Upsilon^2/r_a^2*B1(nt+mt+2)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt) ...
                                        +1/r_a^2*(-k^2*r_a^2-1i*k*r_a*(nt+mt)+mt*nt) ...
                                                *(dradtheta^2*B1(nt+mt)+dradphi^2/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                    +(RtRt + RtR*dradtheta*(-1i*k+mt/r_a)+RRt*dradtheta*(-1i*k+nt/r_a))*B1(nt+mt) ...
                                    +RpRp/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1)) ...
                                    +dradphi/sin(theta)^2*(RpR*(-1i*k+mt/r_a)+RRp*(-1i*k+nt/r_a)) ...
                                                           *(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1));
                        end
                    case 'BGC'
                        if mt+nt == 2
                            temp2 = RR*((1-k^2*Upsilon^2)*B1(2) ...
                                       -Upsilon^2/r_a^2*B1(4)+k^2*Upsilon^2*cos(theta)^2*B1(2) ...
                                       +1/r_a^2*(k^2*r_a^2+1) ...
                                                *(dradtheta^2*B1(2)+dradphi^2/sin(theta)^2*(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3))) ...
                                        -1i*k*r_a) ...
                                   +(RtRt+RtR*dradtheta*(-1i*k+1/r_a)+RRt*dradtheta*(1i*k+1/r_a))*B1(2) ...
                                   +RpRp/sin(theta)^2*(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3)) ...
                                   +dradphi/sin(theta)^2*(RpR*(-1i*k+1/r_a)+RRp*(1i*k+1/r_a)) ...
                                                         *(B2(1)-Upsilon^2/r_a^2*cos(theta)^2*B2(3));
                        else
                            temp2 = RR*(-1i*k*r_a*(nt-mt)*B1(nt+mt-1)+(mt*nt-k^2*Upsilon^2)*B1(nt+mt) ...
                                       +1i*k*Upsilon^2/r_a*(nt-mt)*B1(nt+mt+1)-mt*nt*Upsilon^2/r_a^2*B1(nt+mt+2)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt) ...
                                       +1/r_a^2*(k^2*r_a^2-1i*k*r_a*(nt-mt)+mt*nt) ...
                                                *(dradtheta^2*B1(nt+mt)+dradphi^2/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                   +(RtRt+RtR*dradtheta*(-1i*k+mt/r_a)+RRt*dradtheta*(1i*k+nt/r_a))*B1(nt+mt) ...
                                   +RpRp/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1)) ...
                                   +dradphi/sin(theta)^2*(RpR*(-1i*k+mt/r_a)+RRp*(1i*k+nt/r_a)) ...
                                                         *(B2(nt+mt-1)-Upsilon^2/r_a^2*cos(theta)^2*B2(nt+mt+1));
                        end
                end
                for n = 1:N
                    for m = 1:N
                        if Dt(n,nt)*D(m,mt) ~= 0
                            switch formulation
                                case {'PGU', 'BGU'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*sin(theta)*J_3*J_2*wt*exp(-2*1i*k*r_a);
                                case {'PGC', 'BGC'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*sin(theta)*J_3*J_2*wt;
                            end
                        end
                    end
                end
            end
        end
    end   
    
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    counter = 1;
    for n = 1:N
        for m = 1:N
            indices = counter:(counter+n_en^2-1);
            A_inf_values_temp(indices) = reshape(temp(:,:,n,m),n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctrXiEta+(noSurfDofs*(n-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctrXiEta+(noSurfDofs*(m-1)),n_en,2);
            counter = counter + n_en^2;
        end
    end
    
    A_inf_values(:,e) = A_inf_values_temp;
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noSurfDofs*(N-1);

A_inf = sparse(spIdxRow,spIdxCol,A_inf_values);

% keyboard
if min(size(A_inf)) < noDofs_new
    A_inf(noDofs_new, noDofs_new) = 0;
end

newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));

