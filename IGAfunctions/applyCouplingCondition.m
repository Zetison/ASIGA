function [F, A] = applyCouplingCondition(varCol,alpha_s_arr,beta_s_arr,x_0,k_wn,P_0,omega,rho_f,shift)

if nargin < 9
    shift = 0;
end

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
l = varCol.nurbs.number(3);

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;

weights = varCol.weights;
controlPts = varCol.controlPts;

gluedNodes = varCol.gluedNodes;
noDofs = varCol.noDofs;
noCtrlPts = varCol.noCtrlPts;
noDofs_tot = varCol.noDofs_tot;


% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.


% Check orientation of NURBS object (assuming the object is orientable)
rightHandedOrientation = findOrientation(varCol);

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
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (solidXiEtaMesh == gluedNodes{i}(j));
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
% parfor e = 1:solidNoElemsXiEta
for e = 1:solidNoElemsXiEta
    idXi = solidIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = solidIndexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
    solidSctrXiEta3D = [solidSctrXiEta, solidSctrXiEta+noCtrlPts, solidSctrXiEta+2*noCtrlPts];
    fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:))+noDofs;          %  element scatter vector

    pts = controlPts(solidSctrXiEta,:);
    A1_e = zeros(n_en,3*n_en);
    A2_e = zeros(3*n_en, n_en);
    
    F1_e = zeros(3*n_en,1);
    F2_e = zeros(n_en,1);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(solidNodes));

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        
        if rightHandedOrientation
            normal = crossProd/norm(crossProd);
        else
            normal = -crossProd/norm(crossProd);
        end
        
        v = R_fun*pts;

        A1_e = A1_e + kron(kron(normal',R_fun),R_fun')*norm(crossProd) * J_2 * wt;  
        A2_e = A2_e + kron(R_fun, kron(normal,R_fun'))*norm(crossProd) * J_2 * wt;
        
        k_vec = -[k_wn*cos(beta_s_arr)*sin(alpha_s_arr);
                  k_wn*sin(beta_s_arr)*ones(1,length(alpha_s_arr));
                  k_wn*cos(beta_s_arr)*cos(alpha_s_arr)];
        P_inc = P_0*exp(-1i*dot2(k_vec,x_0)).*exp(1i*dot2(k_vec,v'));
        dP_inc = 1i*dot2(k_vec,normal).*P_inc;
        
        F1_e = F1_e + P_inc*kron(normal,R_fun')*norm(crossProd) * J_2 * wt;
        
        F2_e = F2_e + dP_inc*R_fun'*norm(crossProd) * J_2 * wt;
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

A1 = sparse(spIdxRow1+shift,spIdxCol1+shift,A1values);
A2 = sparse(spIdxRow2+shift,spIdxCol2+shift,A2values);

if min(size(A1)) < noDofs_tot
    A1(noDofs_tot,noDofs_tot) = 0;
end
if min(size(A2)) < noDofs_tot
    A2(noDofs_tot,noDofs_tot) = 0;
end

A = rho_f*omega^2*(A1 + A2);

F = -rho_f*omega^2*vectorAssembly(F1values,indices1+shift,noDofs_tot);
F = F + vectorAssembly(F2values,indices2+shift,noDofs_tot);
