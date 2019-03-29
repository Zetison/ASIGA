% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

if ~exist('F','var')
    F = zeros(noDofs,1);        % external force vector
elseif length(F) ~= noDofs
    F = zeros(noDofs,1);        % external force vector
end

n = n_solid;
m = m_solid;
l = l_solid;

p = p_solid;
q = q_solid;
r = r_solid;

Xi = solid.knots{1};
Eta = solid.knots{2};
Zeta = solid.knots{3};


solidNodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        solidNodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
        counter = counter + 1;
    end
end
fluidNodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        fluidNodes(counter) = n*(j-1) + i;
        counter = counter + 1;
    end
end

[solidXiEtaMesh, solidIndexXiEta, solidNoElemsXiEta, ...
 ~, ~] = generateIGA2DMesh(Xi, Eta, p, q, n, m);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes_solid)
    parentIdx = gluedNodes_solid{i}(1);
    for j = 2:length(gluedNodes_solid{i})
        indices = (solidXiEtaMesh == gluedNodes_solid{i}(j));
        solidXiEtaMesh(indices) = parentIdx;
    end
end

n_en = (p+1)*(q+1);
A1values = zeros(3*n_en^2,solidNoElemsXiEta);
A2values = zeros(3*n_en^2,solidNoElemsXiEta);

spIdxRow1 = zeros(3*n_en^2,solidNoElemsXiEta);
spIdxCol1 = zeros(3*n_en^2,solidNoElemsXiEta);

spIdxRow2 = zeros(3*n_en^2,solidNoElemsXiEta);
spIdxCol2 = zeros(3*n_en^2,solidNoElemsXiEta);

F1values = zeros(3*n_en,solidNoElemsXiEta);
F2values = zeros(n_en,solidNoElemsXiEta);

indices1 = zeros(3*n_en,solidNoElemsXiEta);
indices2 = zeros(n_en,solidNoElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
parfor e = 1:solidNoElemsXiEta
% for e = 1:solidNoElemsXiEta
    idXi = solidIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = solidIndexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
    solidSctrXiEta3D = [solidSctrXiEta, solidSctrXiEta+noCtrlPts_solid, solidSctrXiEta+2*noCtrlPts_solid];
    fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:))+noDofs_solid;          %  element scatter vector

    pts = controlPts_solid(solidSctrXiEta,:);
    A1_e = zeros(n_en,3*n_en);
    A2_e = zeros(3*n_en, n_en);
    
    F1_e = zeros(3*n_en,1);
    F2_e = zeros(n_en,1);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights_solid(solidNodes));

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        
        normal = crossProd/norm(crossProd);

        v = R_fun*pts;
        x = v(1);
        y = v(2);
        z = v(3);

        A1_e = A1_e + kron(kron(normal',R_fun),R_fun')*norm(crossProd) * J_2 * wt;  
        A2_e = A2_e + kron(R_fun, kron(normal,R_fun'))*norm(crossProd) * J_2 * wt;
         
        F1_e = F1_e + P_inc(x,y,z)*kron(normal,R_fun')*norm(crossProd) * J_2 * wt;
        
        F2_e = F2_e + dP_inc(x,y,z, normal)*R_fun'*norm(crossProd) * J_2 * wt;
    end
    spIdxRow1(:,e) = copyVector(fluidSctrXiEta,3*n_en,1);
    spIdxCol1(:,e) = copyVector(solidSctrXiEta3D,n_en,2);
    
    spIdxRow2(:,e) = copyVector(solidSctrXiEta3D,n_en,1);
    spIdxCol2(:,e) = copyVector(fluidSctrXiEta,3*n_en,2);
    
    A1values(:,e) = reshape(A1_e,3*n_en^2,1);
    A2values(:,e) = reshape(A2_e,3*n_en^2,1);
    
    indices1(:,e) = solidSctrXiEta3D';
    indices2(:,e) = fluidSctrXiEta';
    
    F1values(:,e) = F1_e;
    F2values(:,e) = F2_e;
end

A1 = sparse(spIdxRow1,spIdxCol1,A1values);
A2 = sparse(spIdxRow2,spIdxCol2,A2values);

if min(size(A1)) < noDofs_tot
    A1(noDofs_tot,noDofs_tot) = 0;
end
if min(size(A2)) < noDofs_tot
    A2(noDofs_tot,noDofs_tot) = 0;
end

A = A + rho_f*omega^2*(A1 + A2);

F = -rho_f*omega^2*vectorAssembly(F1values,indices1,noDofs_tot);
F = F + vectorAssembly(F2values,indices2,noDofs_tot);

