function nurbs = leastSquares1D(Xi,p_xi,f,d,weights,dofsToRemove,values)
n_xi = length(Xi) - (p_xi+1);
if nargin < 5 || any(isnan(weights))
    weights = ones(1,n_xi);
end
grev = aveknt(Xi, p_xi+1);
controlPts = [grev; weights];
nurbs = createNURBSobject(controlPts,Xi);

varCol.dimension = 3;
varCol.nurbs = nurbs;
varCol = convertNURBS(varCol);
varCol = generateIGAmesh(varCol);
varCol = findDofsToRemove(varCol);

degree = varCol.nurbs{1}.degree;

index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;

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

[Q1D,W1D] = gaussTensorQuad(64);

%% Build global matrices
% for e = 1:noElems
parfor e = 1:noElems

    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = elRange{1}(index(e,1),:);

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);

    xi = parent2ParametricSpace(Xi_e, Q1D);
    I = findKnotSpans(degree, xi(1,:), knots);

    R = NURBSbasis(I,xi,degree,knots,wgts);
    
    J_2 = 0.5*(Xi_e(2)-Xi_e(1));
    f_e = zeros(n_en,d);
    m_e = zeros(n_en);
    f_qps = f(xi);
    
    for gp = 1:size(W1D,1)
        
        m_e = m_e + R{1}(gp,:)'*R{1}(gp,:) * J_2 * W1D(gp);  
        
        f_e = f_e + R{1}(gp,:)'*f_qps(gp,:)* J_2 * W1D(gp);
    end
    spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
    spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
    Mvalues(:,e) = reshape(m_e, sizeMe, 1);
    
    F_indices(:,e) = sctr';
    Fvalues(:,e,:) = f_e;
end

%% Collect data into global matrices (and load vector)
noDofs = noDofs/3;
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

