function F = applyNeumannCondition_CoupledProblem(varCol,omega,rho_f,no_angles,shift)

Xi = varCol.knotVecs{1}{1};
Eta = varCol.knotVecs{1}{2};

p_xi = varCol.patches{1}.nurbs.degree(1);
p_eta = varCol.patches{1}.nurbs.degree(2);

n_xi = varCol.patches{1}.nurbs.number(1);
n_eta = varCol.patches{1}.nurbs.number(2);
n_zeta = varCol.patches{1}.nurbs.number(3);

elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};

weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;

dp_inc = varCol.dp_inc;
p_inc = varCol.p_inc;


gluedNodes = varCol.gluedNodes;
noDofs = varCol.noDofs;
noDofs_tot = varCol.noDofs_tot;

d = varCol.dimension;

% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

solidNodes = zeros(1,n_xi*n_eta);
counter = 1;
for j = 1:n_eta
    for i = 1:n_xi
        solidNodes(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
        counter = counter + 1;
    end
end
fluidNodes = zeros(1,n_xi*n_eta);
counter = 1;
for j = 1:n_eta
    for i = 1:n_xi
        fluidNodes(counter) = n_xi*(j-1) + i;
        counter = counter + 1;
    end
end

[solidXiEtaMesh, solidIndexXiEta, solidNoElemsXiEta] = generateIGA2DMesh(Xi, Eta, p_xi, p_eta, n_xi, n_eta);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (solidXiEtaMesh == gluedNodes{i}(j));
        solidXiEtaMesh(indices) = parentIdx;
    end
end

n_en = (p_xi+1)*(p_eta+1);

F1values = zeros(d*n_en,solidNoElemsXiEta,no_angles);
F2values = zeros(n_en,solidNoElemsXiEta,no_angles);

indices1 = zeros(d*n_en,solidNoElemsXiEta);
indices2 = zeros(n_en,solidNoElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
% warning('parfor is not being used')
% for e = 1:solidNoElemsXiEta
parfor e = 1:solidNoElemsXiEta
    idXi = solidIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = solidIndexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
    solidSctrXiEtadD = zeros(d*length(solidSctrXiEta),1);
    for i = 1:d
        solidSctrXiEtadD(i:d:d*n_en) = d*(solidSctrXiEta-1)+i;
    end
    fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:))+noDofs;          %  element scatter vector

    pts = controlPts(solidSctrXiEta,:);
    
    F1_e = zeros(d*n_en,no_angles);
    F2_e = zeros(n_en,no_angles);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis_old(xi, eta, p_xi, p_eta, Xi, Eta, weights(solidNodes));

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        J_1 = norm(crossProd);
        n = crossProd/J_1;
        
        X = R_fun*pts;
                
        F1_e = F1_e + kron(R_fun',n)*p_inc(X).'*J_1 * J_2 * wt;
        F2_e = F2_e + R_fun'*dp_inc(X,n.').'*J_1 * J_2 * wt;
    end
    
    indices1(:,e) = solidSctrXiEtadD';
    indices2(:,e) = fluidSctrXiEta';
    
    F1values(:,e,:) = F1_e;
    F2values(:,e,:) = F2_e;
end

for alpha_s_Nr = 1:no_angles
    F(:,alpha_s_Nr) = -vectorAssembly(F1values(:,:,alpha_s_Nr),indices1+shift,noDofs_tot);
    F(:,alpha_s_Nr) = F(:,alpha_s_Nr) + 1/(rho_f*omega^2)*vectorAssembly(F2values(:,:,alpha_s_Nr),indices2+shift,noDofs_tot);
end

