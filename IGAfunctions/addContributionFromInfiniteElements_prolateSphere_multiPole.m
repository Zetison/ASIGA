
zeta1Nodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        zeta1Nodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
        counter = counter + 1;
    end
end

[XiEtaMesh, indexXiEta, noElemsXiEta, ...
 ~, ~] = generateIGA2DMesh(Xi, Eta, p, q, n, m);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (XiEtaMesh == gluedNodes{i}(j));
        XiEtaMesh(indices) = parentIdx;
    end
end

temp_n = n;
temp_m = m;
temp_l = l;

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
digits(50)
H1_temp = inv(vpa(S1));
H2_temp = inv(vpa(S2));

H1 = zeros(N);
H2 = zeros(N);
for i = 1:N
    for j = 1:N
        H1(i,j) = H1_temp(i,j);
        H2(i,j) = H2_temp(i,j);
    end
end



% Evaluate analytic integrals in ``radial'' direction. Note that the last
% two integrals (I1(end) and I1(end-1), I2(end) and I2(end-1)) will be redundant for the cases 'BGC'
% and 'BGU'
switch infiniteElementFormulation
    case {'PGU', 'BGU'}
        I1 = zeros(2*N+4,1);
        I2 = zeros(2*N+4,1);
        for n = 1:2*N+4
            I1(n) = 1/R_a^(n-1)*expIntegral(n,-2*1i*k_wn*R_a);
            I2(n) = expIntegral2(n-1,R_a,k_wn,f);
        end
    case {'PGC', 'BGC',}
        I1 = zeros(2*N+4,1);
        I2 = zeros(2*N+4,1);
        for n = 1:2*N+4
            I1(n) = 1/R_a^(n-1)*1/(n-1);
        end
        % Note that I1(1) (= inf) will not be used.
        for n = 1:2*N+4
            I2(n) = fracIntegral2(n-1,R_a,f);
        end
end
% keyboard
% f = evalin(symengine, 'float(int(cos(2*k*r)/r^3, r=R..infinity))');
% matlabFunction(f, 'file','otherFunctions/expIntegral_syms')
% 
% syms rr
% expInt = matlabFunction(cos(2*k_wn*rr)/rr^n);
% % f = exp(2*1i*k_wn*rr)/rr^n;
% int(expInt, rr, R_a, inf)


%% Extend global matrix
A(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;

n_en = (p+1)*(q+1);

spIdxRow = zeros((N*n_en)^2,noElemsXiEta);
spIdxCol = zeros((N*n_en)^2,noElemsXiEta);
A_inf_values = zeros((N*n_en)^2,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
% parfor e = 1:noElemsXiEta
for e = 1:noElemsXiEta
    idXi = indexXiEta(e,1);   % the indexXiEta matrix is made in generateIGA2DMesh
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    
    J1_AB = zeros(n_en);
    J2_AB = zeros(n_en);
    J3_AB = zeros(n_en);
    J4_AB = zeros(n_en);
    J5_AB = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        
%         [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
%         J = pts'*[dRdxi' dRdeta' dRdzeta'];        

        [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));
        J = pts'*[dRdxi' dRdeta'];
        
        v = evaluateNURBS(nurbs,[xi eta 1]);
        x = v(1);
        y = v(2);
        z = v(3);
        
        [radius, theta, phi] = evaluateProlateCoords(x,y,z,f);
        
%         DXDP = dXdP(radius,theta,phi,f);
        DPDX = dPdX(x,y,z,f,radius,theta);
        
%         dRdP = DXDP(2:3,:)*dRdX;
        J3 = DPDX(2:3,:)*J(:,1:2);
        J_3 = abs(det(J3));
        
        dRdtheta = dRdeta*(J3(1,2))^(-1);
        dRdphi = dRdxi*(J3(2,1))^(-1);
                
        J1_AB = J1_AB  + R_fun'*R_fun*sin(theta)                       *J_3*J_2*wt;  
        J2_AB = J2_AB  + dRdtheta'*dRdtheta*sin(theta)                 *J_3*J_2*wt;  
        J3_AB = J3_AB  + dRdphi'*dRdphi/sin(theta)                     *J_3*J_2*wt;  
        J4_AB = J4_AB  + dRdphi'*dRdphi*cos(theta)^2/sin(theta)        *J_3*J_2*wt;  
        J5_AB = J5_AB  + R_fun'*R_fun*cos(theta)^2*sin(theta)          *J_3*J_2*wt;  
%         
    end   
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    for j = 1:N
        for i = 1:N
            switch infiniteElementFormulation
                case {'PGC', 'PGU'}
                    n = i + 2;
                case {'BGC', 'BGU'}
                    n = i;
            end
            m = j;
            temp = zeros(n_en);
            for im = 1:N
                for jm = 1:N
                    switch infiniteElementFormulation
                        case {'PGC', 'PGU'}
                            nm = im + 2;
                        case {'BGC', 'BGU'}
                            nm = im;
                    end
                    mm = jm;
                    switch infiniteElementFormulation
                        case {'PGU', 'BGU'}
                            if mm == 1 && nm == 1
                                blf  = J1_AB*(-2*1i*k_wn*I1(1) + (1+f^2*k_wn^2)*I1(2) + 2*1i*k_wn*f^2*I1(3) - f^2*I1(4) - 1i*k_wn*exp(2*1i*k_wn*R_a)) ...
                                      + J2_AB*I1(2) + J3_AB*I2(1) - J4_AB*f^2*I2(3) + J5_AB*k_wn^2*f^2*I1(2);
                            else
                                blf  = J1_AB*(-2*k_wn^2*I1(nm+mm-2) - 1i*k_wn*(nm+mm)*I1(nm+mm-1) + (nm*mm + f^2*k_wn^2)*I1(nm+mm) ...
                                        + 1i*k_wn*f^2*(nm+mm)*I1(nm+mm+1) - nm*mm*f^2*I1(nm+mm+2)) + J2_AB*I1(nm+mm) ...
                                        + J3_AB*I2(nm+mm-1) - J4_AB*f^2*I2(nm+mm+1) + J5_AB*k_wn^2*f^2*I1(nm+mm);
                            end
                        case {'PGC', 'BGC'}
                            if mm == 1 && nm == 1
                                slf  = J1_AB*((1-f^2*k_wn^2)*I2(1) - f^2*I1(4) - 1i*k_wn) ...
                                       + J2_AB*I1(2) + J3_AB*I2(1) - J4_AB*f^2*I2(3) + J5_AB*k_wn^2*f^2*I1(2);
                            else
                                slf  = J1_AB*(1i*k_wn*(mm-nm)*I1(nm+mm-1) + (nm*mm - f^2*k_wn^2)*I1(nm+mm) ...
                                      - 1i*k_wn*f^2*(mm-nm)*I1(nm+mm+1) - nm*mm*f^2*I1(nm+mm+2)) + J2_AB*I1(nm+mm) ...
                                      + J3_AB*I2(nm+mm-1) - J4_AB*f^2*I2(nm+mm+1) + J5_AB*k_wn^2*f^2*I1(nm+mm);
                            end
                    end
                    switch infiniteElementFormulation
                        case {'PGU', 'BGU'}
                            temp = temp + H1(i,im)*exp(-1i*k_wn*rm(i))*H2(j,jm)*exp(-1i*k_wn*rm(j))*blf;
                        case {'PGC', 'BGC'}
                            temp = temp + H1(i,im)*exp(-1i*k_wn*rm(i))*H2(j,jm)*exp(1i*k_wn*rm(j))*slf;
                    end
                end
            end
            indices = (1:n_en^2)+(n_en^2*N*(j-1) + n_en^2*(i-1));
            A_inf_values_temp(indices) = reshape(temp,n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctrXiEta+(noDofs*(i-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctrXiEta+(noDofs*(j-1)),n_en,2);
        end
    end
    A_inf_values(:,e) = A_inf_values_temp;
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noDofs*(N-1);
if ~memoryCheap
    A_inf = sparse(spIdxRow,spIdxCol,A_inf_values);
    
    F(noDofs_new) = 0; % To get correct dimension of F
    if min(size(A_inf)) < noDofs_new
        A_inf(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;
    end
    A = A + A_inf;    
end



newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));

dofsToRemove = sort(unique([dofsToRemove newDofsToRemove]));


n = temp_n;
m = temp_m;
l = temp_l;


