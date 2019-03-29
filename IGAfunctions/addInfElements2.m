function [A_inf, newDofsToRemove] = addInfElements2(varCol, k, Upsilon, r_a)

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
l = varCol.nurbs.number(3);
gluedNodes = varCol.gluedNodes;
N = varCol.N;
infiniteElementFormulation = varCol.infiniteElementFormulation;
% rm = varCol.rm;

noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;

zeta1Nodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        zeta1Nodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
        counter = counter + 1;
    end
end

[XiEtaMesh, indexXiEta, noElemsXiEta] = generateIGA2DMesh(Xi, Eta, p, q, n, m);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (XiEtaMesh == gluedNodes{i}(j));
        XiEtaMesh(indices) = parentIdx;
    end
end

%% Find coefficients for radial shape functions in infinite elements
% S1 = zeros(N);
% S2 = zeros(N);
% for j = 1:N
%     for jm = 1:N
%         switch infiniteElementFormulation
%             case {'PGC', 'PGU'}
%                 n = j+2;
%             case {'BGC', 'BGU'}
%                 n = j;
%         end
%         m = j;
%         S1(j,jm) = rm(jm)^(-n);
%         S2(j,jm) = rm(jm)^(-m);
%     end
% end
% H1_temp = inv(S1);
% H2_temp = inv(S2);
% 
% 
% H1 = zeros(N);
% H2 = zeros(N);
% for i = 1:N
%     for j = 1:N
%         H1(i,j) = H1_temp(i,j);
%         H2(i,j) = H2_temp(i,j);
%     end
% end
H1 = varCol.D;
H2 = varCol.Dt;


%% Evaluate analytic integrals in ``radial'' direction. 
% Note that the last two integrals (I1(end) and I1(end-1),
% I2(end) and I2(end-1)) will be redundant for the cases 'BGC' and 'BGU'
B1 = zeros(2*N+4,1);
B2 = zeros(2*N+3,1);
for n = 1:2*N+4
    B1(n) = radialIntegral2(n, r_a, k, Upsilon, infiniteElementFormulation, 1);
    if n < 2*N+4
        B2(n) = radialIntegral2(n, r_a, k, Upsilon, infiniteElementFormulation, 2);
    end
end

%% Calculate contribution from infinite elements
n_en = (p+1)*(q+1);

spIdxRow = zeros((N*n_en)^2,noElemsXiEta);
spIdxCol = zeros((N*n_en)^2,noElemsXiEta);
A_inf_values = zeros((N*n_en)^2,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
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
        
        [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));
        J = pts'*[dRdxi' dRdeta'];
        
        v = R_fun*pts;
        x = v(1);
        y = v(2);
        z = v(3);
        
        [r, theta, ~, T] = evaluateProlateCoords(x,y,z,Upsilon);
        chi = r;
        DPDX = dPdX(x,y,z,Upsilon,r,theta,T);
        
        J3 = DPDX(2:3,:)*J(:,1:2);
        J_3 = abs(det(J3));
        
        dRdP = J3'\[dRdxi; dRdeta];
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);
%         dRdtheta = dRdeta*(J3(1,2))^(-1);
%         dRdphi = dRdxi*(J3(2,1))^(-1);

        J1_AB = J1_AB  + R_fun'*R_fun*sin(theta)                       *J_3*J_2*wt;  
        J2_AB = J2_AB  + dRdtheta'*dRdtheta*sin(theta)                 *J_3*J_2*wt;  
        J3_AB = J3_AB  + dRdphi'*dRdphi/sin(theta)                     *J_3*J_2*wt;  
        J4_AB = J4_AB  + dRdphi'*dRdphi*cos(theta)^2/sin(theta)        *J_3*J_2*wt;  
        J5_AB = J5_AB  + R_fun'*R_fun*cos(theta)^2*sin(theta)          *J_3*J_2*wt;  
    end   
    
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    for m = 1:N
        for n = 1:N
            temp = zeros(n_en);
            for nt = 1:N
                for mt = 1:N
%                     blf  = J1_AB*(-2*k^2*chi^2*I1(nt+mt) - 1i*k*chi*(nt+mt+2)*I1(nt+mt+1) + ((nt+2)*mt + Upsilon^2*k^2)*I1(nt+mt+2) ...
%                             + 1i*k*Upsilon^2/chi*(nt+mt+2)*I1(nt+mt+3) - (nt+2)*mt*Upsilon^2/chi^2*I1(nt+mt+4)) + J2_AB*I1(nt+mt+2) ...
%                             + J3_AB*I2(nt+mt+1) - J4_AB*Upsilon^2/chi^2*I2(nt+mt+3) + J5_AB*k^2*Upsilon^2*I1(nt+mt+2);
%                             
                    switch infiniteElementFormulation
                        case 'BGU'
                            if mt+nt == 2
                                blf  = J1_AB*(-2*1i*k*chi*B1(1) + (1+Upsilon^2*k^2)*B1(2) ...
                                      + 2*1i*k*Upsilon^2/chi*B1(3) - Upsilon^2/chi^2*B1(4) - 1i*k*exp(2*1i*k*chi)*chi) ...
                                      + J2_AB*B1(2) + J3_AB*B2(1) - J4_AB*Upsilon^2/chi^2*B2(3) + J5_AB*k^2*Upsilon^2*B1(2);
                            else
                                blf  = J1_AB*(-2*k^2*chi^2*B1(nt+mt-2) - 1i*k*chi*(nt+mt)*B1(nt+mt-1) + (nt*mt + Upsilon^2*k^2)*B1(nt+mt) ...
                                        + 1i*k*Upsilon^2/chi*(nt+mt)*B1(nt+mt+1) - nt*mt*Upsilon^2/chi^2*B1(nt+mt+2)) + J2_AB*B1(nt+mt) ...
                                        + J3_AB*B2(nt+mt-1) - J4_AB*Upsilon^2/chi^2*B2(nt+mt+1) + J5_AB*k^2*Upsilon^2*B1(nt+mt);
                            end
                        case 'PGU'
                            blf  = J1_AB*(-2*k^2*chi^2*B1(nt+mt) - 1i*k*chi*(nt+mt+2)*B1(nt+mt+1) + ((nt+2)*mt + Upsilon^2*k^2)*B1(nt+mt+2) ...
                                    + 1i*k*Upsilon^2/chi*(nt+mt+2)*B1(nt+mt+3) - (nt+2)*mt*Upsilon^2/chi^2*B1(nt+mt+4)) + J2_AB*B1(nt+mt+2) ...
                                    + J3_AB*B2(nt+mt+1) - J4_AB*Upsilon^2/chi^2*B2(nt+mt+3) + J5_AB*k^2*Upsilon^2*B1(nt+mt+2);
                        case 'BGC'
                            if mt+nt == 2
                                blf  = J1_AB*((1-Upsilon^2*k^2)*B1(2) - Upsilon^2/chi^2*B1(4) - 1i*k*chi) ...
                                       + J2_AB*B1(2) + J3_AB*B2(1) - J4_AB*Upsilon^2/chi^2*B2(3) + J5_AB*k^2*Upsilon^2*B1(2);
                            else
                                blf  = J1_AB*(1i*k*chi*(mt-nt)*B1(nt+mt-1) + (nt*mt - Upsilon^2*k^2)*B1(nt+mt) ...
                                      - 1i*k*Upsilon^2/chi*(mt-nt)*B1(nt+mt+1) - nt*mt*Upsilon^2/chi^2*B1(nt+mt+2)) + J2_AB*B1(nt+mt) ...
                                      + J3_AB*B2(nt+mt-1) - J4_AB*Upsilon^2/chi^2*B2(nt+mt+1) + J5_AB*k^2*Upsilon^2*B1(nt+mt);
                            end
                        case 'PGC'
                            blf  = J1_AB*(1i*k*chi*(mt-nt-2)*B1(nt+mt+1) + ((nt+2)*mt - Upsilon^2*k^2)*B1(nt+mt+2) ...
                                  - 1i*k*Upsilon^2/chi*(mt-nt-2)*B1(nt+mt+3) - (nt+2)*mt*Upsilon^2/chi^2*B1(nt+mt+4)) + J2_AB*B1(nt+mt+2) ...
                                  + J3_AB*B2(nt+mt+1) - J4_AB*Upsilon^2/chi^2*B2(nt+mt+3) + J5_AB*k^2*Upsilon^2*B1(nt+mt+2);
                    end
                    
                    switch infiniteElementFormulation
                        case {'PGU', 'BGU'}
                            temp = temp + H1(n,nt)*exp(-2*1i*k*chi)*H2(m,mt)*chi*blf;
                        case {'PGC', 'BGC'}
                            temp = temp + H1(n,nt)*H2(m,mt)*chi*blf;
                    end
%                     temp = temp + H1(n,nt)*exp(-2*1i*k*chi)*H2(m,mt)*chi*blf;
                end
            end
            indices = (1:n_en^2)+(n_en^2*N*(m-1) + n_en^2*(n-1));
            A_inf_values_temp(indices) = reshape(temp,n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctrXiEta+(noDofs*(n-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctrXiEta+(noDofs*(m-1)),n_en,2);
        end
    end
    A_inf_values(:,e) = A_inf_values_temp;
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noDofs*(N-1);

A_inf = sparse(spIdxRow,spIdxCol,A_inf_values);

% keyboard
if min(size(A_inf)) < noDofs_new
    A_inf(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;
end

newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));


