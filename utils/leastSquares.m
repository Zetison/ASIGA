function U = leastSquares(varCol,U, type)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.


%% Extract all needed data from options and varCol
switch type
    case 'scalar'
        d = 1;
    case 'gradient'
        d = 3;
    case 'stress'
        d = 6;
end

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);
p_zeta = varCol.nurbs.degree(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;


dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);

sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
F_indices = zeros(n_en,noElems); 
Fvalues   = zeros(n_en,noElems,d); 

[W3D,Q3D] = gaussianQuadNURBS(p_xi+1,p_eta+1,p_zeta+1); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
% for e = 1:noElems
parfor e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    Zeta_e = elRangeZeta(idZeta,:);
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    
    f_e = zeros(n_en,d);
    m_e = zeros(n_en);
    U_sctr = U(sctr);
    
    for gp = 1:size(W3D,1)
        pt = Q3D(gp,:);
        wt = W3D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = parent2ParametricSpace(Zeta_e,pt(3));
        
        [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);
        
        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta; dRdzeta];

        m_e = m_e + R'*R * abs(J_1) * J_2 * wt;  
        
        switch type
            case 'scalar'
                f_gp = R*U_sctr;
                f_e = f_e + R'*f_gp * abs(J_1) * J_2 * wt;
            case 'gradient'
                f_gp = dRdX*U_sctr;
                f_e = f_e + R.'*f_gp.' * abs(J_1) * J_2 * wt;
        end
    end
    spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
    spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
    Mvalues(:,e) = reshape(m_e, sizeMe, 1);
    
    F_indices(:,e) = sctr';
    Fvalues(:,e,:) = f_e;
end

%% Collect data into global matrices (and load vector)
F = vectorAssembly(Fvalues,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);

M(dofsToRemove,:) = [];  
M(:,dofsToRemove) = [];

F(dofsToRemove,:) = [];

u = M\F;

U = zeros(noDofs,d);
U(setdiff(1:noDofs, dofsToRemove'),:) = u;   
varCol.dimension = 1;
U = addSolutionToRemovedNodes_new(U, varCol);


