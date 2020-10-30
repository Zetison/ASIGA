function F = applyHWBC2(varCol,alpha_s_arr,gradient)
error('Depricated, use applyNeumannCondition instead')

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);

gluedNodes = varCol.gluedNodes;

noDofs_tot = varCol.noDofs_tot;

weights = varCol.weights;
controlPts = varCol.controlPts;

rightHandedOrientation = findOrientation(varCol);


% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs_tot,length(alpha_s_arr));        % external force vector


% Apply Neumann boundary condition at surface where zeta = 0
zeta0Nodes = zeros(1,n*m);
counter = 1;
for j = 1:m
    for i = 1:n
        zeta0Nodes(counter) = n*(j-1) + i;
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

n_en = (p+1)*(q+1);
Fvalues = zeros(n_en,noElemsXiEta,length(alpha_s_arr));
indices = zeros(n_en,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
parfor e = 1:noElemsXiEta
% for e = 1:noElemsXiEta
    idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    f_e = zeros(n_en,length(alpha_s_arr));
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2)); 


        v = R_fun*pts;
        
        if rightHandedOrientation
            normal = -crossProd/norm(crossProd);
        else
            normal = crossProd/norm(crossProd);
        end
        
        f_e = f_e + R_fun'*(gradient(v)*normal)*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctrXiEta';
    Fvalues(:,e,:) = f_e;
end

for alpha_s_Nr = 1:length(alpha_s_arr)
    F(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
end
