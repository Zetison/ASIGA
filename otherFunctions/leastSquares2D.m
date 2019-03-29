function nurbs = leastSquares2D(nurbs,f,d)


Xi = nurbs.knots{1};
Eta = nurbs.knots{2};
p_xi = nurbs.degree(1);
p_eta = nurbs.degree(2);
n_xi = nurbs.number(1);
n_eta = nurbs.number(2);
P = nurbs.coeffs;


varCol.dimension = 1;
varCol = convertNURBS(nurbs, varCol);
varCol = generateIGA2DMesh_new(varCol);

index = varCol.patches{1}.index;
noElems = varCol.patches{1}.noElems;
elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};
element = varCol.patches{1}.element;
element2 = element;
weights = varCol.patches{1}.weights;

% dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.patches{1}.noDofs;

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1);

sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
F_indices = zeros(n_en,noElems); 
Fvalues   = zeros(n_en,noElems,d); 

[W2D,Q2D] = gaussianQuadNURBS(p_xi+10,p_eta+10); 
%% Build global matrices
% warning('Not running in Parallel')
% keyboard
% for e = 1:noElems
parfor e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    sctr = element(e,:);
    wgts = weights(element2(e,:));
    
    f_e = zeros(n_en,d);
    m_e = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        
        R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        
        m_e = m_e + R'*R * J_2 * wt;  
        
        f_e = f_e + R'*f(xi,eta)* J_2 * wt;
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
% U = zeros(noDofs,d);

% [xi_remove, dofsToRemove] = getRepeatedKnots(Xi,p_xi);
dofsToRemove = sort([2:n_xi-1, 1:n_xi:n_xi*n_eta, n_xi:n_xi:n_xi*n_eta, n_xi*(n_eta-1)+2:n_xi*n_eta-1]);

F(dofsToRemove,:) = [];
M(dofsToRemove,:) = [];
P = reshape(P,4,n_xi*n_eta);
F = F - M(:,dofsToRemove)*P(1:3,dofsToRemove).';
% U(dofsToRemove,:) = P(1:3,indices).';

M(:,dofsToRemove) = [];

% U(setdiff(1:noDofs,dofsToRemove),:) = M\F;

% controlPts = [U, ones(n_xi,1)].';

% nurbs = createNURBSobject(controlPts,Xi);

P(1:3,setdiff(1:n_xi*n_eta, dofsToRemove')) = (M\F).';
nurbs.coeffs = reshape(P,4,n_xi,n_eta);




