
% Evaluate analytic integrals in ``radial'' direction. Note that the last
% two integrals (E(end) and E(end-1)) will be redundant for the cases 'BGC'
% and 'BGU'
switch infiniteElementFormulation
    case {'PGU', 'BGU'}
        E = zeros(2*N+2,1);
        for nz = 1:2*N+2
            E(nz) = expIntegral(nz,-2*1i*R_a*k_wn);
%             E(nz) = expIntegral(nz,-2*1i*k_wn);
        end
    case {'PG', 'BGC',}
        % The values of the integrals will be explicitly given
end

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

%% Extend global matrix
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
I_AB_values = zeros(n_en^2,noElemsXiEta);
I_ABm_values = zeros(n_en^2,noElemsXiEta);
J_AB_values = zeros(n_en^2,noElemsXiEta);
P_AB_values = zeros(n_en^2,noElemsXiEta);
Q_AB_values = zeros(n_en^2,noElemsXiEta);

spIdxRow = zeros((N*n_en)^2,noElemsXiEta);
spIdxCol = zeros((N*n_en)^2,noElemsXiEta);
A_inf_values = zeros((N*n_en)^2,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
% parfor e = 1:noElemsXiEta
testint = 0;
for e = 1:noElemsXiEta
    idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    
    I_AB_e  = zeros(n_en);
    I_ABm_e = zeros(n_en);
    J_AB_e  = zeros(n_en);
    Q_AB_e  = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

        J = pts'*[dRdxi' dRdeta' R_fun'];
        J_1 = det(J);
        crossProd = cross(J(:,1), J(:,2)); 
%         kappa = norm(crossProd)/abs(J_1);

        barh_xi   = norm(J(:,1));
        barh_eta  = norm(J(:,2));
        barh_zeta = norm(J(:,3));

        e_xi   = J(:,1)/barh_xi;
        e_eta  = J(:,2)/barh_eta;
        e_zeta = J(:,3)/barh_zeta;

        E_xieta   = e_xi.'*e_eta;
        E_xizeta  = e_xi.'*e_zeta;
        E_etazeta = e_eta.'*e_zeta;
        
        Delta = 1 - E_xieta^2 - E_xizeta^2 - E_etazeta^2 + 2*E_xieta*E_xizeta*E_etazeta;

        a = 1/Delta*[1-E_etazeta^2              E_xizeta*E_etazeta-E_xieta  E_xieta*E_etazeta-E_xizeta;
                     E_xizeta*E_etazeta-E_xieta 1-E_xizeta^2                E_xieta*E_xizeta-E_etazeta;
                     E_xieta*E_etazeta-E_xizeta E_xieta*E_xizeta-E_etazeta  1-E_xieta^2];

        nabla_S_R = 1/barh_xi *kron((a(1,1)*e_xi + a(2,1)*e_eta+a(3,1)*e_zeta),dRdxi) ...
                  + 1/barh_eta*kron((a(1,2)*e_xi + a(2,2)*e_eta+a(3,2)*e_zeta),dRdeta);
              
        a_vec = 1/barh_zeta*(a(1,3)*e_xi + a(2,3)*e_eta + a(3,3)*e_zeta); 

        %%%%%%%%
%         v = J(:,3);
%         x = v(1);
%         y = v(2);
%         z = v(3);
%         radius = sqrt(x^2+y^2+z^2);
%         theta = acos(z/radius);
%         phi = atan2(y,x);
%         if phi < 0
%             phi = phi + 2*pi;
%         end
%         testint  = testint  + exp(theta+phi)               *abs(J_1)*J_2*wt; 
%         testint  = testint  + exp(theta+phi)                *norm(crossProd)*J_2*wt;  
        
%         a_vec = a(1,3)*e_xi + a(2,3)*e_eta + a(3,3)*e_zeta; 
% 
%         1 + norm(a_vec)^2 - 2*a_vec.'*cross(e_xi,e_eta)/norm(cross(e_xi,e_eta))*kappa
        
%         (E_xieta^2-1)/Delta
        
        
%         dRdX = J'\[dRdxi; dRdeta; dRdzeta];
        
        I_AB_e  = I_AB_e  + R_fun'*R_fun                 *abs(J_1)*J_2*wt;  
        I_ABm_e = I_ABm_e + R_fun'*R_fun*norm(a_vec)^2   *abs(J_1)*J_2*wt;  
        J_AB_e  = J_AB_e  + nabla_S_R'*nabla_S_R         *abs(J_1)*J_2*wt;  
        Q_AB_e  = Q_AB_e  + R_fun'*(a_vec'*nabla_S_R)    *abs(J_1)*J_2*wt; 
%         
%         I_AB_e  = I_AB_e  + R_fun'*R_fun                 *norm(crossProd)*J_2*wt;  
%         I_ABm_e = I_ABm_e + R_fun'*R_fun*norm(a_vec)^2   *norm(crossProd)*J_2*wt;  
%         J_AB_e  = J_AB_e  + nabla_S_R'*nabla_S_R         *norm(crossProd)*J_2*wt;  
%         Q_AB_e  = Q_AB_e  + R_fun'*(a_vec'*nabla_S_R)    *norm(crossProd)*J_2*wt;  
    end   
    P_AB_e = Q_AB_e'; %% Eller omvendt???????????????????????????????????????????????????
    switch infiniteElementFormulation
        % Create continuous basis functions by scaling
        case {'PGU', 'BGU'}
            I_AB_e  = exp(-2*1i*R_a*k_wn)*I_AB_e;  
            I_ABm_e = exp(-2*1i*R_a*k_wn)*I_ABm_e;  
            J_AB_e  = exp(-2*1i*R_a*k_wn)*J_AB_e;  
            Q_AB_e  = exp(-2*1i*R_a*k_wn)*Q_AB_e;  
            P_AB_e  = exp(-2*1i*R_a*k_wn)*P_AB_e;
    end
    for j = 1:N
        for i = 1:N
            switch infiniteElementFormulation
                case {'PG', 'PGU'}
                    nz = i + 2;
                case {'BGC', 'BGU'}
                    nz = i;
            end
            mz = j;
            switch infiniteElementFormulation
                case {'PGU', 'BGU'}
                    if mz == 1 && nz == 1
                        temp  = J_AB_e*E(2) ...
                              + P_AB_e*(1i*R_a*k_wn*E(1) - E(2)) ...
                              + Q_AB_e*(1i*R_a*k_wn*E(1) - E(2)) ...
                              + I_ABm_e*(-2*1i*R_a*k_wn*E(1) + E(2)) ...
                              - 1i*k_wn/2*(I_ABm_e + I_AB_e)*exp(2*1i*R_a*k_wn);
                    else
                        temp  = J_AB_e*E(mz+nz) ...
                              + P_AB_e*(1i*R_a*k_wn*E(mz+nz-1) - mz*E(mz+nz)) ...
                              + Q_AB_e*(1i*R_a*k_wn*E(mz+nz-1) - nz*E(mz+nz)) ...
                              + I_ABm_e*(-R_a^2*k_wn^2*E(mz+nz-2) - 1i*R_a*k_wn*(nz+mz)*E(mz+nz-1) + nz*mz*E(mz+nz)) ...
                              - I_AB_e*k_wn^2*E(mz+nz-2);
                    end
                case {'PG', 'BGC'}
                    if mz == 1 && nz == 1
                        error('This case is not implemented as it will only work for the unit sphere anyways.')
                    else
                        temp  = J_AB_e/(mz+nz-1) ...
                              + P_AB_e*(1i*R_a*k_wn/(mz+nz-2) - mz/(mz+nz-1)) ...
                              + Q_AB_e*(-1i*R_a*k_wn/(mz+nz-2) - nz/(mz+nz-1)) ...
                              + I_ABm_e*(R_a^2*k_wn^2/(mz+nz-3) - 1i*R_a*k_wn*(nz-mz)/(mz+nz-2) + nz*mz/(mz+nz-1)) ...
                              - I_AB_e*k_wn^2/(mz+nz-3);
                    end
            end
            indices = (1:n_en^2)+(n_en^2*N*(j-1) + n_en^2*(i-1));
            A_inf_values(indices,e) = reshape(temp,n_en^2,1);
            spIdxRow(indices,e) = copyVector(sctrXiEta+(noDofs*(i-1)),n_en,1);
            spIdxCol(indices,e) = copyVector(sctrXiEta+(noDofs*(j-1)),n_en,2);
        end
    end
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





