function [M, F] = bestApproximationSurf(varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.


%% Extract all needed data from options and varCol

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;


dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1);
v_values   = zeros(n_en,noElems,3); 

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 

%% Find nodes at which to evaluate the exact solution
% for e = 1:noElems
parfor e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    
    v_e = zeros(size(W2D,1),3);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        
        R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
        
        v_e(gp,:) = R*pts;
    end
    v_values(:,e,:) = v_e;
end
v_values = reshape(v_values, size(W2D,1)*noElems, 3);
analytic_values = varCol.analytic(v_values);
analytic_values = reshape(analytic_values, size(W2D,1), noElems);


%% Build global matrices
sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
F_values   = zeros(n_en,noElems); 
F_indices = zeros(n_en,noElems); 

parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    
    F_e = zeros(n_en,1);
    m_e = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

        J = pts'*[dRdxi' dRdeta'];
        crossProd = cross(J(:,1),J(:,2));
        J_1 = norm(crossProd);

        m_e = m_e + R'*R * abs(J_1) * J_2 * wt;  
        
        p_analytic = analytic_values(gp,e);
        F_e = F_e + p_analytic*R.' * abs(J_1) * J_2 * wt; 
    end
    spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
    spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
    Mvalues(:,e) = reshape(m_e, sizeMe, 1);
    
    F_indices(:,e) = sctr';
    F_values(:,e) = F_e;
end

%% Collect data into global matrices (and load vector)
F = vectorAssembly(F_values,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);



