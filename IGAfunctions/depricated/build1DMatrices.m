
% initialization
error('Depricated. Use buildMatrices instead')
% The program will execute faster if the K matrix is not made sparse untill
% right before the system is solved. This is due to the fact that indexing
% is slow, new memory must at all time be allocated as the matrix grows.

U = zeros(noDofs,1);        % displacement vector
n_en = p+1;
sizeOfk_e = n_en^2;
indices = zeros(n_en,noElems); 
spIdxRow = zeros(sizeOfk_e,noElems);
spIdxCol = zeros(sizeOfk_e,noElems);
Kvalues = zeros(sizeOfk_e,noElems); 
if buildMassMatrix
    Mvalues = zeros(sizeOfk_e,noElems); 
end
Fvalues = zeros(n_en,noElems); 



% Gauss quadrature rule
[W1D,Q1D] = gaussianQuadNURBS(p+1); 

% Assembling system of equation
% Stiffness matrix and external force vector
% Loop over elements (knot spans)
for e = 1:noElems
    idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
    
    % calculate element ranges
    % elRangeU, elRangeV and elRangeW comes from generateIGA3DMesh
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    
    J_2 = 0.5*(Xi_e(2)-Xi_e(1));
    
    sctr = element(e,:);          %  element scatter vector
    pts = controlPts(sctr);
    
    k_e = zeros(n_en);
    if buildMassMatrix
        m_e = zeros(n_en);
    end
    f_e = zeros(n_en,1);
    
    for gp = 1:size(W1D,1)
        pt = Q1D(gp);
        wt = W1D(gp);

        xi = parent2ParametricSpace(Xi_e,  pt);
        
        [R_fun, dRdxi] = NURBS1DBasis(xi, p, Xi, weights);
        
        J = dRdxi*pts;
        J_1 = abs(J);
        dRdx = dRdxi/J;
        
        k_e = k_e + dRdx' * dRdx * abs(J_1) * J_2 * wt; 
        if buildMassMatrix
            m_e = m_e + R_fun'*R_fun * abs(J_1) * J_2 * wt;      
        end
    end
    spIdxRow(:,e) = copyVector(sctr,n_en,1);
    spIdxCol(:,e) = copyVector(sctr,n_en,2);
    Kvalues(:,e) = reshape(k_e, sizeOfk_e, 1);

    if buildMassMatrix
        Mvalues(:,e) = reshape(m_e, sizeOfk_e, 1);
    end
    indices(:,e) = sctr';
end
K = E*sparse(spIdxRow,spIdxCol,Kvalues);
if buildMassMatrix
    M = sparse(spIdxRow,spIdxCol,Mvalues);
end
if min(size(K)) < noDofs
    K(noDofs,noDofs) = 0;
    M(noDofs,noDofs) = 0;
end

if buildMassMatrix
    M = rho_s*M;
end
