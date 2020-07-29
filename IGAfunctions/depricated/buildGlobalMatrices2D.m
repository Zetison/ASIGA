function [K, M, F] = buildGlobalMatrices2D(varCol, newOptions)
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
buildMassMatrix = options.buildMassMatrix;
applyBodyLoading = options.applyBodyLoading;
operator = options.operator;

Xi = varCol.patches{1}.nurbs.knots{1};
Eta = varCol.patches{1}.nurbs.knots{2};
p = varCol.patches{1}.nurbs.degree(1);
q = varCol.patches{1}.nurbs.degree(2);

index = varCol.patches{1}.index;
noElems = varCol.patches{1}.noElems;
elRangeXi = varCol.patches{1}.elRange{1};
elRangeEta = varCol.patches{1}.elRange{2};
element = varCol.patches{1}.element;
weights = varCol.patches{1}.weights;
controlPts = varCol.patches{1}.controlPts;

noCtrlPts = varCol.patches{1}.noCtrlPts;
noDofs = varCol.patches{1}.noDofs;

if strcmp(options.operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = (p+1)*(q+1);
sizeKe = (d*n_en)^2;

spIdxRow = zeros((d*n_en)^2,noElems);
spIdxCol = zeros((d*n_en)^2,noElems);
Kvalues = zeros((d*n_en)^2,noElems); 

if buildMassMatrix
    Mvalues = zeros(sizeKe,noElems); 
end
if applyBodyLoading
    F_indices = zeros(d*n_en,noElems); 
    Fvalues = zeros(d*n_en,noElems); 
    force = varCol.f;
else
    force = NaN;
end

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
% [W2D,Q2D] = gaussianQuadNURBS(60,60); 

%% Build global matrices
parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element(e,:),:); % New
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    k_e = zeros(d*n_en);
    if buildMassMatrix
        m_e = zeros(n_en);
    end
    if applyBodyLoading
        f_e = zeros(d*n_en,1);
    end
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, wgts);
        
        J = pts'*[dRdxi' dRdeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta];
        
        switch operator
            case 'linearElasticity'
                B = strainDispMatrix2d(n_en,dRdX);
                k_e = k_e + B.' * C * B * abs(J_1) * J_2 * wt; 
            case 'Laplace'
                k_e = k_e + dRdX'*dRdX* abs(J_1) * J_2 * wt;  
        end
        if buildMassMatrix
            m_e = m_e + R'*R * abs(J_1) * J_2 * wt;  
        end
        
        if applyBodyLoading
            v = R*pts;
            f_gp = force(v(1),v(2));
            f_e = f_e + kron(R.', f_gp.') * abs(J_1) * J_2 * wt;
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
    Kvalues(:,e) = reshape(temp, sizeKe, 1);

    if buildMassMatrix
        temp = zeros(d*n_en,d*n_en);
        for i = 1:d
            temp(i:d:end, i:d:end) = m_e;
        end
        Mvalues(:,e) = reshape(temp, sizeKe, 1);
    end
    if applyBodyLoading
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
% 
% k_ee = 0;
% for xi = [-1/sqrt(3),1/sqrt(3)]
%     x = (xi+1)/2;
%     for eta = [-1/sqrt(3),1/sqrt(3)]
%         y = (eta+1)/2;
%         k_ee = k_ee + (x^2+y^2+(x+y)^2/2)*0.25;
%     end
% end
% k_ee
