function p_h = calculateScatteredPressure2(varCol, U, k, P_far, noWavesVec)

n_res = 6;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);
p_zeta = varCol.nurbs.degree(3);
n_xi = varCol.nurbs.number(1);
n_eta = varCol.nurbs.number(2);
n_zeta = varCol.nurbs.number(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;


% Check orientation of NURBS object (assuming the object is orientable)
rightHandedOrientation = findOrientation(varCol);

% Find elements on the inner surface for evaluation of
% backscattered pressure in far field
innerSurfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(2) == 1
        innerSurfaceElements = [innerSurfaceElements e];
    end
end

noElemsXi = length(unique(varCol.nurbs.knots{1})) - 1;
noElemsEta = length(unique(varCol.nurbs.knots{2})) - 1;

noWavesInXiDir = noWavesVec(1);
noWavesInEtaDir = noWavesVec(2);

noGpXi = ceil(n_res*noWavesInXiDir/noElemsXi);
noGpEta = ceil(n_res*noWavesInEtaDir/noElemsEta);

if noGpXi <= p_xi+1
    noGpXi = p_xi+1;
end
if noGpEta <= p_eta+1
    noGpEta = p_eta+1;
end

[W2D,Q2D] = gaussianQuadNURBS(noGpXi,noGpEta);  

if size(P_far,2) < size(U,2)
    error('The number of angles in P_far must be greater or equal to the number of angles in U')
end
if size(U,3) > 1
    error('This is not implemented, but ready to use if changed to matrix3Dprod(A,B)...')
end

Upsilon = varCol.Upsilon;
noDofs = varCol.noDofs;
rm = varCol.rm;
infiniteElementFormulation = varCol.infiniteElementFormulation;

N = varCol.N;
D = varCol.D;
S1 = zeros(N);
S2 = zeros(N);
for j = 1:N
    for jm = 1:N
        switch infiniteElementFormulation
            case {'PGC', 'PGU'}
                n = j+2;
            case {'BGC', 'BGU'}
                n = j;
        end
        m = j;
        S1(j,jm) = rm(jm)^(-n);
        S2(j,jm) = rm(jm)^(-m);
    end
end
H1_temp = inv(S1);
H2_temp = inv(S2);


H1 = zeros(N);
H2 = zeros(N);
for i = 1:N
    for j = 1:N
        H1(i,j) = H1_temp(i,j);
        H2(i,j) = H2_temp(i,j);
    end
end
p_h = zeros(1,size(P_far,2), size(P_far,3));
p_h_gp_arr = zeros(length(innerSurfaceElements)*size(W2D,1),1);
dp_h_gp_arr = zeros(length(innerSurfaceElements)*size(W2D,1),1);
theta_arr = zeros(length(innerSurfaceElements)*size(W2D,1),1);
phi_arr = zeros(length(innerSurfaceElements)*size(W2D,1),1);
% 
% x_hat = zeros(size(P_far));
% for i = 1:size(P_far,2)
%     for j = 1:size(P_far,3)
%         x_hat(:,i,j) = P_far(:,i,j)/norm(P_far(:,i,j));
%     end
% end
% scalingFact = exp(1i*k*norm(P_far(:,1,1)))/norm(P_far(:,1,1));

counter = 1;
for i = 1:length(innerSurfaceElements) %
% parfor i = 1:length(innerSurfaceElements) %[8 7 3 4 5 6 2 1]% 

    e = innerSurfaceElements(i);
    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);

    pts = controlPts(sctr,:);

    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));        
        zeta = 1;

        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);

        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        crossProd = cross(J(:,1),J(:,2));
        if rightHandedOrientation
            normal = crossProd/norm(crossProd);
        else
            normal = -crossProd/norm(crossProd);
        end
        
        X = R_fun*pts;
        diffr = zeros(size(P_far));
        diffr(1,:,:) = X(1)-P_far(1,:,:);
        diffr(2,:,:) = X(2)-P_far(2,:,:);
        diffr(3,:,:) = X(3)-P_far(3,:,:);
        
        norm_rm_r = norm2(diffr);    
        greensFunc = exp(1i*k*norm_rm_r)./(4*pi*norm_rm_r);    
        dgreensFunc = greensFunc.*(1i*k*norm_rm_r - 1)./(norm_rm_r.^2) .*dot3(diffr, normal);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x = X(1);
        y = X(2);
        z = X(3);

        [r, theta, phi, T, c_1, c_2] = evaluateProlateCoords(x,y,z,Upsilon);

        chi = r;

        DPDX = dPdX(x,y,z,Upsilon,r,theta,T);

        J3 = DPDX(2:3,:)*J;

        [dchidtheta, dchidphi] = dchidP(x,y,z,Upsilon,c_1,c_2,J,J3);

        DXDP = dXdP(r,theta,phi,Upsilon);

        h_r = norm(DXDP(:,1));
        h_theta = norm(DXDP(:,2));
        h_phi = norm(DXDP(:,3));

        e_r = DXDP(:,1)/h_r;
        e_theta = DXDP(:,2)/h_theta;
        e_phi = DXDP(:,3)/h_phi;
        p_h_gp = zeros(1,size(U,2));
        dp_h_gp = zeros(1,size(U,2));
        infElVar = varCol.infElVar;
        switch infElVar
            case 1
                for m = 1:N
                    temp = 0;
                    temp2 = 0;
                    for mt = 1:N
                        temp = temp + (1i*k - mt/chi)*H2(m,mt)/chi^mt;
                        temp2 = temp2 + H2(m,mt)/chi^mt;
                    end
                    temp = temp*exp(1i*k*(chi-rm(m)));
                    temp2 = temp2*exp(1i*k*(chi-rm(m)));
                    dp_h_gp = dp_h_gp + temp*R_fun*U(sctr + noDofs*(m-1),:);
                    p_h_gp = p_h_gp + temp2*R_fun*U(sctr + noDofs*(m-1),:);
                end
            case 2
                for m = 1:N
                    temp = 0;
                    temp2 = 0;
                    for mt = 1:N
                        temp = temp + (1i*k - mt/chi)*D(m,mt);
                        temp2 = temp2 + D(m,mt);
                    end
                    dp_h_gp = dp_h_gp + temp*R_fun*U(sctr + noDofs*(m-1),:);
                    p_h_gp = p_h_gp + temp2*R_fun*U(sctr + noDofs*(m-1),:);
                end
            case 3
                for m = 1:N
                    temp = 0;
                    temp2 = 0;
                    for mt = 1:N
                        temp = temp + (1i*k - mt/chi)*D(m,mt);
                        temp2 = temp2 + D(m,mt);
                    end
                    dp_h_gp = dp_h_gp + temp*R_fun*U(sctr + noDofs/n_zeta*(m-1),:);
                    p_h_gp = p_h_gp + temp2*R_fun*U(sctr + noDofs/n_zeta*(m-1),:);
                end
        end
        dp_h_gp = dp_h_gp*(dot(normal,e_r)/h_r - dot(normal,e_theta)/h_theta*dchidtheta ...
                                               - dot(normal,e_phi)/h_phi*dchidphi);
        p_h_gp_arr(counter) = p_h_gp;
        dp_h_gp_arr(counter) = dp_h_gp;
        theta_arr(counter) = theta;
        phi_arr(counter) = phi;
        counter = counter + 1;
%         dp_h_gp = dp_h_gp*(-0.020839827711645 + 0.197893864481245i)/(0.085526631791482 + 0.080731020892433i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_h = p_h + (p_h_gp.*dgreensFunc - dp_h_gp.*greensFunc)* norm(crossProd) * J_2 * wt; 
%         keyboard
%         x_hat_Dot_ny = dot3(x_hat, normal);
%         x_hat_Dot_y = dot3(x_hat, X');
%         p_h = p_h + scalingFact*1/(2*pi)*(dp_h_gp + 1i*k*p_h_gp*x_hat_Dot_ny).*exp(-1i*k*x_hat_Dot_y)*norm(crossProd) * J_2 * wt;  
%         keyboard
    end
%     if length(P_far) > 3
%         disp(['Completed ' num2str(i) ' out of ' num2str(length(innerSurfaceElements)) ' Elapsed time: ' num2str(toc)])
%     end
end
% if false
%     indices = find(abs(phi_arr - phi_arr(1)) < 1e-13);
% %     plot(theta_arr(indices), real(dp_h_gp_arr(indices)))
%     plot(theta_arr(indices), real(p_h_gp_arr(indices)))
% else
%     [phi_arr, I] = sort(phi_arr);
%     theta_arr = theta_arr(I);
%     p_h_gp_arr = p_h_gp_arr(I);
%     dp_h_gp_arr = dp_h_gp_arr(I);
% 
% 
%     Phi = reshape(phi_arr, noGpEta*noElemsEta, noGpXi*noElemsXi);
%     Theta = reshape(theta_arr, noGpEta*noElemsEta, noGpXi*noElemsXi);
% %     Z = reshape(p_h_gp_arr, noGpEta*noElemsEta, noGpXi*noElemsXi);
%     Z = reshape(dp_h_gp_arr, noGpEta*noElemsEta, noGpXi*noElemsXi);
%     for i = 1:noGpXi*noElemsXi
%         [Theta(:,i), I] = sort(Theta(:,i));
%         Phi(:,i) = Phi(I,i);
%         Z(:,i) = Z(I,i);
%     end
%     P_0 = varCol.P_0;
%     R_o = varCol.R_o;
%     dpExact = P_0*(1i*k-1/R_o)*exp(1i*k*(r-R_o))/R_o;
%     surf(Phi,Theta,real(Z));
% %     surf(Phi,Theta,real(Z-dpExact));
%     xlabel('Phi')
%     ylabel('Theta')
% end
% 
% 
% keyboard





