function [A_inf, newDofsToRemove] = infElementsNonSepGeom2(varCol, k, Upsilon,chimin,chimax,scatterer)

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;
noElems = varCol.noElems;
index = varCol.index;
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

N = varCol.N;
infiniteElementFormulation = varCol.infiniteElementFormulation;

noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

%Use chebychev polynomials
D = generateCoeffMatrix(N);
Dt = D;
useApproxRadialIntegrals = 2;
B1splines = {};
B2splines = {};
B1cheb = {};
B2cheb = {};

if useApproxRadialIntegrals == 1
    p_splineApprox = 2;
    B1splines = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, infiniteElementFormulation, 1),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
    B2splines = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, infiniteElementFormulation, 2),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
elseif useApproxRadialIntegrals == 2
    B1cheb = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                infiniteElementFormulation, 1), 1e-6);
    end
    B2cheb = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                infiniteElementFormulation, 2), 1e-6);
    end
end
disp('done calculating integrals')
%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

spIdxRow = zeros((N*n_en)^2,noElems);
spIdxCol = zeros((N*n_en)^2,noElems);
A_inf_values = zeros((N*n_en)^2,noElems);
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1);

% maxchi_recorded = -inf;
% minchi_recorded = inf;

% for e = 1:noElems
parfor e = 1:noElems
%     e
%     noElems
    idXi = index(e,1); 
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);
    
    n_en = length(sctr);
    pts = controlPts(sctr,:);
    
    temp = zeros(n_en,n_en,N,N);
    for gp = 1:size(W2D,1)
%         tic
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
        
        J = pts'*[dRdxi' dRdeta'];
        
        X = R_fun*pts;
        x = X(1);
        y = X(2);
        z = X(3);
%         hold on
%         plot3(v(1),v(2),v(3),'*')
        [r, theta, phi, T, c_1, c_2] = evaluateProlateCoords(x,y,z,Upsilon);
        objFun = @(prm) norm(evaluateNURBS(scatterer,[prm(2), prm(3)]) - [sqrt(prm(1)^2-Upsilon^2)*sin(theta)*cos(phi);
                                                                          sqrt(prm(1)^2-Upsilon^2)*sin(theta)*sin(phi);
                                                                          prm(1)*cos(theta)]);
        output = fminsearchbnd(objFun, ([chimin, 0, 0]+[chimax, 1, 1])/2, [chimin, 0, 0],[chimax, 1, 1]);
        r_R = output(1);
        X = [sqrt(r_R^2-Upsilon^2)*sin(theta)*cos(phi);
             sqrt(r_R^2-Upsilon^2)*sin(theta)*sin(phi);
             r_R*cos(theta)];
        x = X(1);
        y = X(2);
        z = X(3);
        chi = r_R;
%         if abs(chi-r) > 1e-13
%             keyboard
%         end
        if chi > chimax || chi < chimin
            error('chi is not in range [chimin chimax]')
        end
%         if chi < minchi_recorded
%             minchi_recorded = chi;
%         end
%         if chi > maxchi_recorded
%             maxchi_recorded = chi;
%         end
        DPDX = dPdX(x,y,z,Upsilon,r,theta,T);
        
        J3 = DPDX(2:3,:)*J;
        J_3 = abs(det(J3));
        
        [dchidtheta, dchidphi] = dchidP(x,y,z,Upsilon,c_1,c_2,J,J3);
        
        dRdP = J3'\[dRdxi; dRdeta];
        
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);
%         toc
%         tic

        if useApproxRadialIntegrals == 1
            B1 = zeros(2*N+4,1);
            xi_pt = (chi-chimin)/(chimax-chimin);
            for i = 1:2*N+4
                B1(i) = evaluateNURBS(B1splines{i}, xi_pt);
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                B2(i) = evaluateNURBS(B2splines{i}, xi_pt);
            end
        elseif useApproxRadialIntegrals == 2
            B1 = zeros(2*N+4,1);
            x_cheb = (2*chi-chimax-chimin)/(chimax-chimin);
            
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
                B1(i) = radialIntegral2(i, chi, k, Upsilon, infiniteElementFormulation, 1);
            end
            B2 = zeros(2*N+3,1);
            for i = 1:2*N+3
                B2(i) = radialIntegral2(i, chi, k, Upsilon, infiniteElementFormulation, 2);
            end      
        end
%         keyboard
        RR = R_fun'*R_fun;  
        RtRt = dRdtheta'*dRdtheta;  
        RpRp = dRdphi'*dRdphi;  
        
        RtR = dRdtheta'*R_fun;  
        RRt = R_fun'*dRdtheta; 
        RpR = dRdphi'*R_fun; 
        RRp = R_fun'*dRdphi;
        
        for nt = 1:N
            for mt = 1:N
                switch infiniteElementFormulation
                    case 'PGU'
                        temp2 = RR*(-2*k^2*chi^2*B1(nt+mt)-1i*k*chi*(nt+mt+2)*B1(nt+mt+1)+(mt*(nt+2)+k^2*Upsilon^2)*B1(nt+mt+2) ...
                                   +1i*k*Upsilon^2/chi*(nt+mt+2)*B1(nt+mt+3)-mt*(nt+2)*Upsilon^2/chi^2*B1(nt+mt+4)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt+2) ...
                                   +1/chi^2*(-k^2*chi^2-1i*k*chi*(nt+mt+2)+mt*(nt+2)) ...
                                           *(dchidtheta^2*B1(nt+mt+2)+dchidphi^2/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3)))) ...
                               +(RtRt+RtR*dchidtheta*(-1i*k+mt/chi)+RRt*dchidtheta*(-1i*k+(nt+2)/chi))*B1(nt+mt+2) ...
                               +RpRp/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3)) ...
                               +dchidphi/sin(theta)^2*(RpR*(-1i*k+mt/chi)+RRp*(-1i*k+(nt+2)/chi)) ...
                                                     *(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3));
                                                 
%                         temp2 = RR*(-2*k^2*chi^2*B1(nt+mt)-1i*k*chi*(nt+mt+2)*B1(nt+mt+1)+(mt*(nt+2)+k^2*Upsilon^2)*B1(nt+mt+2) ...
%                                    +1i*k*Upsilon^2/chi*(nt+mt+2)*B1(nt+mt+3)-mt*(nt+2)*Upsilon^2/chi^2*B1(nt+mt+4)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt+2)) ...
%                                +RtRt*B1(nt+mt+2) ...
%                                +RpRp/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3));
                                                 
%                         temp2 = RR*(-k^2*chi^2*B1(nt+mt)-1i*k*chi*(nt+mt+2)*B1(nt+mt+1) ...
%                                     +mt*(nt+2)*B1(nt+mt+2)+k^2*Upsilon^2*B1(nt+mt+2) ...
%                                     +1i*k*Upsilon^2/chi*(nt+mt+2)*B1(nt+mt+3) - mt*(nt+2)*Upsilon^2/chi^2*B1(nt+mt+4)) ...
%                                 +RtRt*B1(nt+mt+2)+RtR*dchidtheta*(-1i*k*B1(nt+mt+2)+mt/chi*B1(nt+mt+2)) ...
%                                 +RRt*dchidtheta*(-1i*k*B1(nt+mt+2)+(nt+2)/chi*B1(nt+mt+2)) ...
%                                 +RR*dchidtheta^2*(-k^2*B1(nt+mt+2)-1i*k*(nt+mt+2)/chi*B1(nt+mt+2) ...
%                                                             +1/chi^2*mt*(nt+2)*B1(nt+mt+2)) ...
%                                 +1/sin(theta)^2*(RpRp*B2(nt+mt+1) ...
%                                                  +RpR*dchidphi*(-1i*k*B2(nt+mt+1)+mt/chi*B2(nt+mt+1)) ...
%                                                  +RpR*dchidphi*(-1i*k*B2(nt+mt+1)+(nt+2)/chi*B2(nt+mt+1)) ...
%                                                  +RR*dchidphi^2*(-k^2*B2(nt+mt+1)-1i*k*(mt+nt+2)/chi*B2(nt+mt+1) ...
%                                                                  +mt*(nt+2)/chi^2*B2(nt+mt+1))) ...
%                                 -Upsilon^2*cos(theta)^2/(chi^2*sin(theta)^2)*(RpRp*B2(nt+mt+3) ...
%                                                  +RpR*dchidphi*(-1i*k*B2(nt+mt+3)+mt/chi*B2(nt+mt+3)) ...
%                                                  +RpR*dchidphi*(-1i*k*B2(nt+mt+3)+(nt+2)/chi*B2(nt+mt+3)) ...
%                                                  +RR*dchidphi^2*(-k^2*B2(nt+mt+3)-1i*k*(mt+nt+2)/chi*B2(nt+mt+3) ...
%                                                                  +mt*(nt+2)/chi^2*B2(nt+mt+3)))...
%                                 -k^2*RR*(chi^2*B1(nt+mt)-Upsilon^2*cos(theta)^2*B1(nt+mt+2));

%                         temp2 = RR*(-k^2*(chi^2*B1(nt+mt)-Upsilon^2*B1(nt+mt+2)) ...
%                                     -(mt+nt+2)/chi*(chi^2*B1(nt+mt+1)-Upsilon^2*B1(nt+mt+3))*1i*k ...
%                                     +mt*(nt+2)/chi^2*(chi^2*B1(nt+mt+2)-Upsilon^2*B1(nt+mt+4))) ...
%                                 +RtRt*B1(nt+mt+2)+(RtR*(-1i*k+mt/chi)+RRt*(-1i*k+(nt+2)/chi))*dchidtheta*B1(nt+mt+2) ...
%                                 +RR*(-k^2-(mt+nt+2)/chi*1i*k+mt*(nt+2)/chi^2)*dchidtheta^2*B1(nt+mt+2) ...
%                                 +1/chi^2*1/sin(theta)^2*( RpRp+(RpR*(-1i*k+mt/chi)+RRp*(-1i*k+(nt+2)/chi))*dchidphi ...
%                                                          +RR*(-k^2-(mt+nt+2)/chi*1i*k+mt*(nt+2)/chi^2)*dchidphi^2)...
%                                                         *(chi^2*B2(nt+mt+1)-Upsilon^2*cos(theta)^2*B2(nt+mt+3)) ...
%                                 -k^2*(chi^2*B1(nt+mt)-Upsilon^2*cos(theta)^2*B1(nt+mt+2))*RR;
                            
                    case 'PGC'
                        temp2 = RR*(-1i*k*chi*(nt-mt+2)*B1(nt+mt+1)+(mt*(nt+2)-k^2*Upsilon^2)*B1(nt+mt+2) ...
                                   +1i*k*Upsilon^2/chi*(nt-mt+2)*B1(nt+mt+3)-mt*(nt+2)*Upsilon^2/chi^2*B1(nt+mt+4)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt+2) ...
                                   +1/chi^2*(k^2*chi^2-1i*k*chi*(nt-mt+2)+mt*(nt+2)) ...
                                           *(dchidtheta^2*B1(nt+mt+2)+dchidphi^2/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3)))) ...
                               +(RtRt+RtR*dchidtheta*(-1i*k+mt/chi)+RRt*dchidtheta*(1i*k+(nt+2)/chi))*B1(nt+mt+2) ...
                               +RpRp/sin(theta)^2*(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3)) ...
                               +dchidphi/sin(theta)^2*(RpR*(-1i*k+mt/chi)+RRp*(1i*k+(nt+2)/chi)) ...
                                                     *(B2(nt+mt+1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+3));
                    case 'BGU'
                        if mt+nt == 2
                            temp2 = RR*(-2*1i*k*chi*B1(1)+(1+k^2*Upsilon^2)*B1(2) ...
                                        +2*1i*k*Upsilon^2/chi*B1(3)-Upsilon^2/chi^2*B1(4)+k^2*Upsilon^2*cos(theta)^2*B1(2) ...
                                        +1/chi^2*(-k^2*chi^2-2*1i*k*chi+1) ...
                                           *(dchidtheta^2*B1(2)+dchidphi^2/sin(theta)^2*(B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3))) ...
                                        -1i*k*exp(2*1i*k*chi)*chi) ...
                                    +(RtRt + RtR*dchidtheta*(-1i*k+1/chi)+RRt*dchidtheta*(-1i*k+1/chi))*B1(2) ...
                                    +RpRp/sin(theta)^2*(B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3)) ...
                                    +dchidphi/sin(theta)^2*(RpR*(-1i*k+1/chi)+RRp*(-1i*k+1/chi)) ...
                                                           *B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3);
                        else
                            temp2 = RR*(-2*k^2*chi^2*B1(nt+mt-2)-1i*k*chi*(nt+mt)*B1(nt+mt-1)+(mt*nt+k^2*Upsilon^2)*B1(nt+mt) ...
                                        +1i*k*Upsilon^2/chi*(nt+mt)*B1(nt+mt+1)-mt*nt*Upsilon^2/chi^2*B1(nt+mt+2)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt) ...
                                        +1/chi^2*(-k^2*chi^2-1i*k*chi*(nt+mt)+mt*nt) ...
                                                *(dchidtheta^2*B1(nt+mt)+dchidphi^2/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                    +(RtRt + RtR*dchidtheta*(-1i*k+mt/chi)+RRt*dchidtheta*(-1i*k+nt/chi))*B1(nt+mt) ...
                                    +RpRp/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1)) ...
                                    +dchidphi/sin(theta)^2*(RpR*(-1i*k+mt/chi)+RRp*(-1i*k+nt/chi)) ...
                                                           *B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1);
                        end
                    case 'BGC'
                        if mt+nt == 2
                            temp2 = RR*((1-k^2*Upsilon^2)*B1(2) ...
                                       -Upsilon^2/chi^2*B1(4)+k^2*Upsilon^2*cos(theta)^2*B1(2) ...
                                       +1/chi^2*(k^2*chi^2+1) ...
                                                *(dchidtheta^2*B1(2)+dchidphi^2/sin(theta)^2*(B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3))) ...
                                        -1i*k*chi) ...
                                   +(RtRt+RtR*dchidtheta*(-1i*k+1/chi)+RRt*dchidtheta*(1i*k+1/chi))*B1(2) ...
                                   +RpRp/sin(theta)^2*(B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3)) ...
                                   +dchidphi/sin(theta)^2*(RpR*(-1i*k+1/chi)+RRp*(1i*k+1/chi)) ...
                                                         *(B2(1)-Upsilon^2/chi^2*cos(theta)^2*B2(3));
                        else
                            temp2 = RR*(-1i*k*chi*(nt-mt)*B1(nt+mt-1)+(mt*nt-k^2*Upsilon^2)*B1(nt+mt) ...
                                       +1i*k*Upsilon^2/chi*(nt-mt)*B1(nt+mt+1)-mt*nt*Upsilon^2/chi^2*B1(nt+mt+2)+k^2*Upsilon^2*cos(theta)^2*B1(nt+mt) ...
                                       +1/chi^2*(k^2*chi^2-1i*k*chi*(nt-mt)+mt*nt) ...
                                                *(dchidtheta^2*B1(nt+mt)+dchidphi^2/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1)))) ...
                                   +(RtRt+RtR*dchidtheta*(-1i*k+mt/chi)+RRt*dchidtheta*(1i*k+nt/chi))*B1(nt+mt) ...
                                   +RpRp/sin(theta)^2*(B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1)) ...
                                   +dchidphi/sin(theta)^2*(RpR*(-1i*k+mt/chi)+RRp*(1i*k+nt/chi)) ...
                                                         *(B2(nt+mt-1)-Upsilon^2/chi^2*cos(theta)^2*B2(nt+mt+1));
                        end
                end
                for n = 1:N
                    for m = 1:N
                        if D(n,nt)*D(m,mt) ~= 0
                            switch infiniteElementFormulation
                                case {'PGU', 'BGU'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*chi*exp(-2*1i*k*chi)*sin(theta)*J_3*J_2*wt;
                                case {'PGC', 'BGC'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*chi*sin(theta)*J_3*J_2*wt;
                            end
                        end
                    end
                end
            end
        end
    end  
%     keyboard
    
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
    
    A_inf_values(:,e) = A_inf_values_temp;
    
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end
% maxchi_recorded
% minchi_recorded
noDofs_new = noDofs + noDofs*(N-1);

A_inf = sparse(spIdxRow,spIdxCol,A_inf_values);


if min(size(A_inf)) < noDofs_new
    A_inf(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;
end
newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));

