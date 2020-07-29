function [K, M, F] = buildGlobalMatrices1D(varCol, newOptions)
error('Depricated. Use buildMatrices instead')
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.

%% Interpret input arguments

% set default values
options = struct('operator','Laplace',...
                 'fieldDimension',1,...
                 'buildMassMatrix',0,...
                 'applyBodyLoading',0);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(newOptions);
if round(nArgs/2) ~= nArgs/2
	error('Must have propertyName/propertyValue pairs')
end

for pair = reshape(newOptions,2,[]) %# pair is {propName;propValue}
    inpName = pair{1}; %# make case insensitive

    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

%% Extract all needed data from options and varCol
d = options.fieldDimension;

Xi = varCol.patches{1}.nurbs.knots;
p_xi = varCol.patches{1}.nurbs.degree;

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
element = varCol.element;
weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;
f = varCol.f;
noCtrlPts = varCol.patches{1}.noCtrlPts;
noDofs = varCol.patches{1}.noDofs;

if strcmp(options.operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = p_xi+1;

spIdxRow = zeros((d*n_en)^2,noElems);
spIdxCol = zeros((d*n_en)^2,noElems);
Kvalues = zeros((d*n_en)^2,noElems); 

if options.buildMassMatrix
    Mvalues = zeros((d*n_en)^2,noElems); 
end
if options.applyBodyLoading
    F_indices = zeros(d*n_en,noElems); 
    Fvalues = zeros(d*n_en,noElems); 
end

[W1D,Q1D] = gaussianQuadNURBS(p_xi+1); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
for e = 1:noElems
% parfor e = 1:noElems
    idXi = index(e);
    
    Xi_e = elRangeXi(idXi,:);
    
    J_2 = 0.5*(Xi_e(2)-Xi_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr);
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    k_e = zeros(d*n_en);
    if options.buildMassMatrix
        m_e = zeros(n_en);
    end
    if options.applyBodyLoading
        f_e = zeros(d*n_en,1);
    end
    
    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi = parent2ParametricSpace(Xi_e,  pt(1));
        
        [R_fun, dRdxi] = NURBS1DBasis(xi, p_xi, Xi, weights);
        
        J = dRdxi*pts;
        J_1 = J;
        dRdX = dRdxi/J_1;
        
        switch options.operator
            case 'linearElasticity'
                k_e = k_e + C * abs(J_1) * J_2 * wt; 
            case 'laplace'
                k_e = k_e + dRdX'*dRdX* abs(J_1) * J_2 * wt;  
        end
        if options.buildMassMatrix
            m_e = m_e + R_fun'*R_fun * abs(J_1) * J_2 * wt;  
        end
        
        if options.applyBodyLoading
            X = R_fun*pts;
            f_gp = f(X);
            f_e = f_e + kron(R_fun', f_gp) * abs(J_1) * J_2 * wt;
        end
    end
    spIdxRow(:,e) = copyVector(sctr_k_e,d*n_en,1);
    spIdxCol(:,e) = copyVector(sctr_k_e,d*n_en,2);
    temp = zeros(d*n_en,d*n_en);
    for i = 1:d
        for j = 1:d
            temp(i:d:end, j:d:end) = k_e(1+(i-1)*n_en:i*n_en, 1+(j-1)*n_en:j*n_en);
        end
    end
    Kvalues(:,e) = reshape(temp, (d*n_en)^2, 1);
    
    if options.buildMassMatrix
        temp = zeros(d*n_en,d*n_en);
        for i = 1:d
            temp(i:d:end, i:d:end) = m_e;
        end
        Mvalues(:,e) = reshape(temp, (d*n_en)^2, 1);
    end
    if options.applyBodyLoading
        F_indices(:,e) = sctr_k_e';
        Fvalues(:,e) = f_e;
    end
end

%% Collect data into global matrices (and load vector)
if options.applyBodyLoading
    F = vectorAssembly(Fvalues,F_indices,noDofs);
end

K = sparse(spIdxRow,spIdxCol,Kvalues);
if options.buildMassMatrix
    M = sparse(spIdxRow,spIdxCol,Mvalues);
else
    M = [];
end

if min(size(K)) < noDofs
    K(noDofs,noDofs) = 0;
    if options.buildMassMatrix
        M(noDofs,noDofs) = 0;
    end
end
