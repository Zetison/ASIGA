function nurbs = leastSquares1D(Xi,p_xi,f,d,weights,dofsToRemove,values)
n_xi = length(Xi) - (p_xi+1);
if nargin < 5 || any(isnan(weights))
    weights = ones(1,n_xi);
end
if isa(weights,'function_handle')
    error('wrong use of function')
end
grev = aveknt(Xi, p_xi+1);
controlPts = [grev; weights];
nurbs = createNURBSobject(controlPts,Xi);

varCol.dimension = 1;
varCol = convertNURBS(nurbs, varCol);
varCol = generateIGA1DMesh(varCol);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
element = varCol.element;
weights = varCol.patches{1}.weights;

% dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.patches{1}.noDofs;

%% Preallocation and initiallizations
n_en = p_xi+1;

sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
F_indices = zeros(n_en,noElems); 
Fvalues   = zeros(n_en,noElems,d); 

[W1D,Q1D] = gaussianQuadNURBS(64); 

%% Build global matrices
% for e = 1:noElems
parfor e = 1:noElems
    idXi = index(e,1);
    
    Xi_e = elRangeXi(idXi,:);
    
    J_2 = 0.5*(Xi_e(2)-Xi_e(1));
    
    sctr = element(e,:);
    wgts = weights(sctr);
    f_e = zeros(n_en,d);
    m_e = zeros(n_en);
    
    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        
        R = NURBS1DBasis(xi,p_xi,Xi, wgts);
        
        m_e = m_e + R'*R * J_2 * wt;  
        
        f_e = f_e + R'*f(xi)* J_2 * wt;
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
U = zeros(noDofs,d);

if nargin < 6
    [xi_remove, dofsToRemove] = getRepeatedKnots(Xi,p_xi);

    F(dofsToRemove,:) = [];
    M(dofsToRemove,:) = [];
    F = F - M(:,dofsToRemove)*f(xi_remove.');
    M(:,dofsToRemove) = [];
    U(setdiff(1:noDofs,dofsToRemove),:) = M\F;
    U(dofsToRemove,:) = f(xi_remove.');
else
    for i = 1:d
        Mtemp = M;
        dofsToRemoveTemp = dofsToRemove{i};
        Ftemp = F(:,i);
        Ftemp(dofsToRemoveTemp,:) = [];
        Mtemp(dofsToRemoveTemp,:) = [];
        Ftemp = Ftemp - Mtemp(:,dofsToRemoveTemp)*values{i};
        Mtemp(:,dofsToRemoveTemp) = [];
        U(setdiff(1:noDofs,dofsToRemoveTemp),i) = Mtemp\Ftemp;
        U(dofsToRemoveTemp,i) = values{i};
    end
end
    



controlPts = [U, weights].';

nurbs = createNURBSobject(controlPts,Xi);

