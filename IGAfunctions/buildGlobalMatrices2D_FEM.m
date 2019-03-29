function [K, M, F] = buildGlobalMatrices2D_FEM(varCol, newOptions)
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
if nargin > 1
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
end

%% Extract all needed data from options and varCol
d = options.fieldDimension;

Xi = varCol.patches{1}.nurbs.knots{1};
Eta = varCol.patches{1}.nurbs.knots{2};
p_xi = varCol.patches{1}.nurbs.degree(1);
p_eta = varCol.patches{1}.nurbs.degree(2);

index = varCol.patches{1}.index;
noElems = varCol.patches{1}.noElems;
elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};
element = varCol.patches{1}.element;
weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;

j_lag_xi = createIdxMap(p_xi+1);
i_lag_xi = repmat((1:p_xi+1)',1,p_xi);
j_lag_eta = createIdxMap(p_eta+1);
i_lag_eta = repmat((1:p_eta+1)',1,p_eta);
% GLL_p = GLLpoints(p+1);  
% GLL_q = GLLpoints(q+1);  
GLL_p = linspace(-1,1,p_xi+1);  
GLL_q = linspace(-1,1,p_eta+1);  

noCtrlPts = varCol.patches{1}.noCtrlPts;
noDofs = varCol.patches{1}.noDofs;

if strcmp(options.operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1);

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

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
% [W2D,Q2D] = gaussianQuadNURBS(60,60); 

%% Build global matrices
parfor e = 1:noElems
% for e = 1:noElems
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(1+(i-1)*n_en:i*n_en) = sctr+(i-1)*noCtrlPts;
    end
    k_e = zeros(d*n_en);
    if options.buildMassMatrix
        m_e = zeros(d*n_en);
    end
    if options.applyBodyLoading
        f_e = zeros(d*n_en,1);
    end
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);
        
        [L_xi,dL_xi] = lagrangeBasis(pt(1),i_lag_xi,GLL_p,j_lag_xi);
        [L_eta,dL_eta] = lagrangeBasis(pt(2),i_lag_eta,GLL_q,j_lag_eta);
        
        R = kron(L_xi,L_eta);
        dRdxi = kron(dL_xi,L_eta);
        dRdeta = kron(L_xi,dL_eta);
        J = pts'*[dRdxi' dRdeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta];
        
        switch options.operator
            case 'linearElasticity'
                B = strainDispMatrix3d(n_en,dRdX);
                k_e = k_e + B' * C * B * abs(J_1) * wt; 
                if options.buildMassMatrix
                    m_e = m_e + blkdiag(R'*R, R'*R) * abs(J_1) * wt;      
                end
            case 'Laplace'
                k_e = k_e + dRdX'*dRdX* abs(J_1) * wt;  
                if options.buildMassMatrix
                    m_e = m_e + R'*R * abs(J_1) * wt;  
                end
        end
        
        
        if options.applyBodyLoading
            v = R*pts;
            f_gp = varCol.f(v(1),v(2));
            f_e = f_e + kron(f_gp, R') * abs(J_1) * wt;
        end
    end

    spIdxRow(:,e) = copyVector(sctr_k_e,d*n_en,1);
    spIdxCol(:,e) = copyVector(sctr_k_e,d*n_en,2);
    Kvalues(:,e) = reshape(k_e, (d*n_en)^2, 1);

    if options.buildMassMatrix
        Mvalues(:,e) = reshape(m_e, (d*n_en)^2, 1);
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
