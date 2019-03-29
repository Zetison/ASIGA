
zeta1Nodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        zeta1Nodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
        counter = counter + 1;
    end
end

temp_n = n;
temp_m = m;
temp_l = l;

% Evaluate analytic integrals in ``radial'' direction. Note that the last
% two integrals (E(end) and E(end-1)) will be redundant for the cases 'BGC'
% and 'BGU'
switch infiniteElementFormulation
    case {'PGU', 'BGU'}
        I1 = zeros(2*N+4,1);
        I2 = zeros(2*N+4,1);
        for i = 1:2*N+4
            I1(i) = 1/R_a^(i-1)*expIntegral(i,-2*1i*k_wn*R_a);
            I2(i) = expIntegral2(i-1,R_a,k_wn,f);
        end
    case {'PGC', 'BGC',}
        I1 = zeros(2*N+4,1);
        I2 = zeros(2*N+4,1);
        for i = 1:2*N+4
            I1(i) = 1/R_a^(i-1)*1/(i-1);
        end
        % Note that I1(1) (= inf) will not be used.
        for i = 1:2*N+4
            I2(i) = fracIntegral2(i-1,R_a,f);
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

%% Extend global matrix
% A_temp = A;
% 
% A = spalloc(noDofs+(N-1)*noDofs,noDofs+(N-1)*noDofs,nnz(A)+2*N^2*nnz(A(zeta1Nodes,1:noDofs)) + (N-1)^2*nnz(A(zeta1Nodes,zeta1Nodes)));
% 
% A(1:noDofs,1:noDofs) = A_temp;
A(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;

for j = 1:N
    for i = 1:N
        A(zeta1Nodes+(j-1)*noDofs,1:noDofs) = A(zeta1Nodes,1:noDofs);
        A(1:noDofs,zeta1Nodes+(i-1)*noDofs) = A(1:noDofs,zeta1Nodes);
    end
end
for j = 2:N
    for i = 2:N
        A(zeta1Nodes+(j-1)*noDofs,zeta1Nodes+(i-1)*noDofs) = A(zeta1Nodes,zeta1Nodes);
    end
end

n_en = (p+1)*(q+1);

spIdxRow = zeros((N*n_en)^2,noElemsXiEta);
spIdxCol = zeros((N*n_en)^2,noElemsXiEta);
A_inf_values = zeros((N*n_en)^2,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
parfor e = 1:noElemsXiEta
% testint = 0;
% for e = 1:noElemsXiEta
    idXi = indexXiEta(e,1);   % the indexXiEta matrix is made in generateIGA2DMesh
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    
    J1_AB_e = zeros(n_en);
    J2_AB_e = zeros(n_en);
    J3_AB_e = zeros(n_en);
    J4_AB_e = zeros(n_en);
    J5_AB_e = zeros(n_en);
    
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
                
        J1_AB_e = J1_AB_e  + R_fun'*R_fun*sin(theta)                       *J_3*J_2*wt;  
        J2_AB_e = J2_AB_e  + dRdtheta'*dRdtheta*sin(theta)                 *J_3*J_2*wt;  
        J3_AB_e = J3_AB_e  + dRdphi'*dRdphi/sin(theta)                     *J_3*J_2*wt;  
        J4_AB_e = J4_AB_e  + dRdphi'*dRdphi*cos(theta)^2/sin(theta)        *J_3*J_2*wt;  
        J5_AB_e = J5_AB_e  + R_fun'*R_fun*cos(theta)^2*sin(theta)          *J_3*J_2*wt;  
        
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
            switch infiniteElementFormulation
                % Create continuous basis functions by scaling
                case {'PGU', 'BGU'}
                    J1_AB = R_a^(m+n)*exp(-2*1i*k_wn*R_a)*J1_AB_e;  
                    J2_AB = R_a^(m+n)*exp(-2*1i*k_wn*R_a)*J2_AB_e;  
                    J3_AB = R_a^(m+n)*exp(-2*1i*k_wn*R_a)*J3_AB_e;  
                    J4_AB = R_a^(m+n)*exp(-2*1i*k_wn*R_a)*J4_AB_e;  
                    J5_AB = R_a^(m+n)*exp(-2*1i*k_wn*R_a)*J5_AB_e;  
                    if m == 1 && n == 1
                        temp  = J1_AB*(-2*1i*k_wn*I1(1) + (1+f^2*k_wn^2)*I1(2) + 2*1i*k_wn*f^2*I1(3) - f^2*I1(4) - 1i*k_wn*exp(2*1i*k_wn*R_a)) ...
                                + J2_AB*I1(2) + J3_AB*I2(1) - J4_AB*f^2*I2(3) + J5_AB*k_wn^2*f^2*I1(2);
                    else
                        temp  = J1_AB*(-2*k_wn^2*I1(n+m-2) - 1i*k_wn*(n+m)*I1(n+m-1) + (n*m + f^2*k_wn^2)*I1(n+m) ...
                                + 1i*k_wn*f^2*(n+m)*I1(n+m+1) - n*m*f^2*I1(n+m+2)) + J2_AB*I1(n+m) ...
                                + J3_AB*I2(n+m-1) - J4_AB*f^2*I2(n+m+1) + J5_AB*k_wn^2*f^2*I1(n+m);
                    end
                case {'PGC', 'BGC'}
                    % Create continuous basis functions by scaling
                    J1_AB = R_a^(m+n)*J1_AB_e;  
                    J2_AB = R_a^(m+n)*J2_AB_e;  
                    J3_AB = R_a^(m+n)*J3_AB_e;  
                    J4_AB = R_a^(m+n)*J4_AB_e;  
                    J5_AB = R_a^(m+n)*J5_AB_e; 
                    if m == 1 && n == 1
                        temp  = J1_AB*((1-f^2*k_wn^2)*I2(1) - f^2*I1(4) - 1i*k_wn) + J2_AB*I1(2) + J3_AB*I2(1) - J4_AB*f^2*I2(3) + J5_AB*k_wn^2*f^2*I1(2);
                    else
                        temp  = J1_AB*(1i*k_wn*(m-n)*I1(n+m-1) + (n*m - f^2*k_wn^2)*I1(n+m) ...
                              - 1i*k_wn*f^2*(m-n)*I1(n+m+1) - n*m*f^2*I1(n+m+2)) + J2_AB*I1(n+m) ...
                              + J3_AB*I2(n+m-1) - J4_AB*f^2*I2(n+m+1) + J5_AB*k_wn^2*f^2*I1(n+m);
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


