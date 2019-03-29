function A = applyCouplingCondition_new(varCol,shift)

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

gluedNodes = varCol.gluedNodes;
noDofs = varCol.noDofs;
noCtrlPts = varCol.noCtrlPts;
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
A1values = zeros(d*n_en^2,solidNoElemsXiEta);
A2values = zeros(d*n_en^2,solidNoElemsXiEta);

spIdxRow1 = zeros(d*n_en^2,solidNoElemsXiEta);
spIdxCol1 = zeros(d*n_en^2,solidNoElemsXiEta);

spIdxRow2 = zeros(d*n_en^2,solidNoElemsXiEta);
spIdxCol2 = zeros(d*n_en^2,solidNoElemsXiEta);


[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
parfor e = 1:solidNoElemsXiEta
% for e = 1:solidNoElemsXiEta
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
    A1_e = zeros(n_en,d*n_en);
    A2_e = zeros(d*n_en, n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights(solidNodes));

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
        
        normal = crossProd/norm(crossProd);
        
        A1_e = A1_e + kron(kron(normal',R_fun),R_fun')*norm(crossProd) * J_2 * wt;  
        A2_e = A2_e + kron(R_fun, kron(normal,R_fun'))*norm(crossProd) * J_2 * wt;        
    end
    spIdxRow1(:,e) = copyVector(fluidSctrXiEta,d*n_en,1);
    spIdxCol1(:,e) = copyVector(solidSctrXiEtadD,n_en,2);
    
    spIdxRow2(:,e) = copyVector(solidSctrXiEtadD,n_en,1);
    spIdxCol2(:,e) = copyVector(fluidSctrXiEta,d*n_en,2);
    
    temp = zeros(n_en,d*n_en);
    for j = 1:d
        temp(:, j:d:end) = A1_e(:, 1+(j-1)*n_en:j*n_en);
    end
    A1values(:,e) = reshape(temp,d*n_en^2,1);
    
    temp = zeros(d*n_en,n_en);
    for i = 1:d
        temp(i:d:end, :) = A2_e(1+(i-1)*n_en:i*n_en, :);
    end
    A2values(:,e) = reshape(temp,d*n_en^2,1);
end

spIdxRow1 = reshape(spIdxRow1,numel(spIdxRow1),1);
spIdxCol1 = reshape(spIdxCol1,numel(spIdxCol1),1);
A1values = reshape(A1values,numel(A1values),1);

[spIdx1,~,IuniqueIdx1] = unique([spIdxRow1, spIdxCol1]+shift,'rows');
A1values = accumarray(IuniqueIdx1,A1values);


spIdxRow2 = reshape(spIdxRow2,numel(spIdxRow2),1);
spIdxCol2 = reshape(spIdxCol2,numel(spIdxCol2),1);
A2values = reshape(A2values,numel(A2values),1);

[spIdx2,~,IuniqueIdx2] = unique([spIdxRow2, spIdxCol2]+shift,'rows');
A2values = accumarray(IuniqueIdx2,A2values);


A =   sparse(spIdx1(:,1),spIdx1(:,2),A1values,noDofs_tot,noDofs_tot,numel(IuniqueIdx1)) ...
    + sparse(spIdx2(:,1),spIdx2(:,2),A2values,noDofs_tot,noDofs_tot,numel(IuniqueIdx2));
