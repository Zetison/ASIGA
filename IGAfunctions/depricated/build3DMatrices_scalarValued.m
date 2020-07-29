
% initialization
error('Depricated. Use buildMatrices instead')
% The program will execute faster if the K matrix is not made sparse untill
% right before the system is solved. This is due to the fact that indexing
% is slow, new memory must at all time be allocated as the matrix grows.

% U = zeros(noDofs,1);        % displacement vector

n_en = (p+1)*(q+1)*(r+1);
if memoryCheap
    noNonZeroEntries = findNoNonZeroEntries(noDofs,element,noElems,1);
    K = spalloc(noDofs,noDofs,noDofs^2-noNonZeroEntries);
    M = spalloc(noDofs,noDofs,noDofs^2-noNonZeroEntries);
else
    Kvalues = zeros(n_en^2,noElems); 
    Mvalues = zeros(n_en^2,noElems); 
    spIdxRow = zeros(n_en^2,noElems);
    spIdxCol = zeros(n_en^2,noElems);
end


[W3D,Q3D] = gaussianQuadNURBS(p+1,q+1,r+1); 

parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
        
    sctr = element(e,:);
    
    pts = controlPts(sctr,:);
    
    k_e = zeros(n_en);
    f_e = zeros(n_en,1);
    m_e = zeros(n_en);
    
    for gp = 1:size(W3D,1)
        pt = Q3D(gp,:);
        wt = W3D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = parent2ParametricSpace(Zeta_e,pt(3));

        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);

        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta; dRdzeta];
        
        k_e = k_e + dRdX'*dRdX* abs(J_1) * J_2 * wt;  
        m_e = m_e + R_fun'*R_fun * abs(J_1) * J_2 * wt;  
        
    end
    if memoryCheap
        for i = 1:n_en
            ii = sctr(i);
            for j = 1:n_en
                jj = sctr(j);
%                 K(ii,jj) = K(ii,jj) + k_e(i,j);  
%                 M(ii,jj) = M(ii,jj) + m_e(i,j);  
            end
        end 
    else
        Kvalues(:,e) = reshape(k_e, n_en^2, 1);
        Mvalues(:,e) = reshape(m_e, n_en^2, 1);
        spIdxRow(:,e) = copyVector(sctr,n_en,1);
        spIdxCol(:,e) = copyVector(sctr,n_en,2);
    end
end
if ~memoryCheap    
    K = sparse(spIdxRow,spIdxCol,Kvalues);
    M = sparse(spIdxRow,spIdxCol,Mvalues);
    if min(size(K)) < noDofs
        K(noDofs,noDofs) = 0;
        M(noDofs,noDofs) = 0;
    end    
end
M = -k_wn^2*M;
A = K+M;