
% initialization
% The program will execute faster if the K matrix is not made sparse untill
% right before the system is solved. This is due to the fact that indexing
% is slow, new memory must at all time be allocated as the matrix grows.

U = zeros(noDofs,1);        % displacement vector
n_en = (p+1)*(q+1)*(r+1);
if applyBodyLoading
    F = zeros(noDofs,1);        % external force vector
end

if memoryCheap
    noNonZeroEntries = findNoNonZeroEntries(noDofs,element,noElems,3);
    K = spalloc(noDofs,noDofs,noDofs^2-noNonZeroEntries);
    if buildMassMatrix
        M = spalloc(noDofs,noDofs,noDofs^2-noNonZeroEntries);
    end
else
    sizeOfk_e = (3*n_en)^2;
    indices = zeros(3*n_en,noElems); 
    spIdxRow = zeros(sizeOfk_e,noElems);
    spIdxCol = zeros(sizeOfk_e,noElems);
    Kvalues = zeros(sizeOfk_e,noElems); 
    if buildMassMatrix
        Mvalues = zeros(sizeOfk_e,noElems); 
    end
    Fvalues = zeros(3*n_en,noElems); 
end



% Gauss quadrature rule
[W3D,Q3D] = gaussianQuadNURBS(p+1,q+1,r+1); 

% Assembling system of equation
% Stiffness matrix and external force vector
% Loop over elements (knot spans)
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = index(e,2);
    idZeta = index(e,3);
    
    % calculate element ranges
    % elRangeU, elRangeV and elRangeW comes from generateIGA3DMesh
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
    
    sctr = element(e,:);          %  element scatter vector
    sctrB = [sctr sctr+noCtrlPts sctr+2*noCtrlPts]; % scatters a B matrix
    pts = controlPts(sctr,:);
    B = zeros(6,3*n_en);
    k_e = zeros(3*n_en);
    if buildMassMatrix
        m_e = zeros(3*n_en);
    end
    if applyBodyLoading
        f_e = zeros(3*n_en,1);
    end
    
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

        B = strainDispMatrix3d(n_en,dRdX);
        
        k_e = k_e + B' * C * B * abs(J_1) * J_2 * wt; 
        if buildMassMatrix
            m_e = m_e + blkdiag(R_fun'*R_fun, R_fun'*R_fun, R_fun'*R_fun) * abs(J_1) * J_2 * wt;      
        end
        
        if applyBodyLoading
            v = R_fun*pts;
            f_gp = f(v(1),v(2),v(3));
            f_e = f_e + [f_gp(1)*R_fun'; f_gp(2)*R_fun'; f_gp(3)*R_fun'] * abs(J_1) * J_2 * wt;
        end
    end

    if memoryCheap
        for i = 1:3*n_en
            ii = sctrB(i);
            for j = 1:3*n_en
                jj = sctrB(j);
                
%                 M(ii,jj) = M(ii,jj) + m_e(i,j);  
%                 K(ii,jj) = K(ii,jj) + k_e(i,j);
            end
        end
    else
        spIdxRow(:,e) = copyVector(sctrB,3*n_en,1);
        spIdxCol(:,e) = copyVector(sctrB,3*n_en,2);
        Kvalues(:,e) = reshape(k_e, sizeOfk_e, 1);
        
        if buildMassMatrix
            Mvalues(:,e) = reshape(m_e, sizeOfk_e, 1);
        end
        if applyBodyLoading
            indices(:,e) = sctrB';
        end
    end
    if applyBodyLoading
        Fvalues(:,e) = f_e;
    end
end
if ~memoryCheap
    if applyBodyLoading
        F = F + vectorAssembly(Fvalues,indices,noDofs);
    end
    
    K = sparse(spIdxRow,spIdxCol,Kvalues);
    if buildMassMatrix
        M = sparse(spIdxRow,spIdxCol,Mvalues);
    end
    if min(size(K)) < noDofs
        K(noDofs,noDofs) = 0;
        M(noDofs,noDofs) = 0;
    end
end

if buildMassMatrix
    M = rho_s*M;
end
