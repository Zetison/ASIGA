function [K, M, F] = buildGlobalMatrices(varCol, newOptions)
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
buildMassMatrix = options.buildMassMatrix;
applyBodyLoading = options.applyBodyLoading;
operator = options.operator;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
elRangeZeta = varCol.elRange{3};
pIndex = varCol.pIndex;
noDofs = varCol.noDofs;

if strcmp(operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);
sizeKe = (d*n_en)^2;
spIdxRow = zeros(sizeKe,noElems,'uint32');
spIdxCol = zeros(sizeKe,noElems,'uint32');
Kvalues = zeros(sizeKe,noElems); 

if buildMassMatrix
    Mvalues = zeros(sizeKe,noElems); 
end
if applyBodyLoading
    F_indices = zeros(d*n_en,noElems); 
    Fvalues = zeros(d*n_en,noElems); 
    force = varCol.force;
else
    force = NaN;
end

[W3D,Q3D] = gaussianQuadNURBS(p_xi+1,p_eta+1,p_zeta+1); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    Zeta = knotVecs{patch}{3}; % New
    
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    Zeta_e = elRangeZeta(idZeta,:);
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:); % New
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
    
    for gp = 1:size(W3D,1)
        pt = Q3D(gp,:);
        wt = W3D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = parent2ParametricSpace(Zeta_e,pt(3));
        
        [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts); % New
        
        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta; dRdzeta];
        
        switch operator
            case 'linearElasticity'
                B = strainDispMatrix3d(n_en,dRdX);
                k_e = k_e + B.' * C * B * abs(J_1) * J_2 * wt; 
            case 'Laplace'
                k_e = k_e + dRdX'*dRdX* abs(J_1) * J_2 * wt;  
        end
        if buildMassMatrix
            m_e = m_e + R'*R * abs(J_1) * J_2 * wt;  
        end
        
        if applyBodyLoading
            v = R*pts;
            f_gp = force(v(1),v(2),v(3));
            f_e = f_e + kron(R', f_gp) * abs(J_1) * J_2 * wt;
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
if applyBodyLoading
    F = vectorAssembly(Fvalues,F_indices,noDofs);
end

if 0
    primes_noElems = factor(noElems);
    prt = prod(primes_noElems(1:3));
    Kvalues = reshape(Kvalues, sizeKe*prt,noElems/prt);
    Mvalues = reshape(Mvalues, sizeKe*prt,noElems/prt);
    spIdxRow = reshape(spIdxRow, sizeKe*prt,noElems/prt);
    spIdxCol = reshape(spIdxCol, sizeKe*prt,noElems/prt);
    for i = 1:noElems/prt
        [spIdx,~,I] = unique([spIdxRow(:,i), spIdxCol(:,i)],'rows','stable');
        KvaluesTemp = accumarray(I,Kvalues(:,i));
        if buildMassMatrix
            MvaluesTemp = accumarray(I,Mvalues(:,i));
        end
        m = numel(KvaluesTemp);
        if m < sizeKe*prt
            Kvalues(:,i) = [KvaluesTemp; zeros(sizeKe*prt-m,1)];
            spIdxRow(:,i) = [spIdx(:,1); zeros(sizeKe*prt-m,1)];
            spIdxCol(:,i) = [spIdx(:,2); zeros(sizeKe*prt-m,1)];
            if buildMassMatrix
                Mvalues(:,i) = [MvaluesTemp; zeros(sizeKe*prt-m,1)];
            end
        end
    end
%     if noElems > 260000
%         keyboard
%     end
    nnzEntries = sum(~(~spIdxRow(:)));
    spIdxRowUnique = zeros(nnzEntries,1,'uint32');
    spIdxColUnique = zeros(nnzEntries,1,'uint32');
    KvaluesUnique = zeros(nnzEntries,1);
    if buildMassMatrix
        MvaluesUnique = zeros(nnzEntries,1);
    end
    counter = 1;
    for i = 1:noElems/prt
        indices = find(spIdxRow(:,i));
        m = numel(indices);
        spIdxRowUnique(counter:counter+m-1) = spIdxRow(indices,i);
        spIdxColUnique(counter:counter+m-1) = spIdxCol(indices,i);
        KvaluesUnique(counter:counter+m-1) = Kvalues(indices,i);
        if buildMassMatrix
            MvaluesUnique(counter:counter+m-1) = Mvalues(indices,i);
        end
        counter = counter + m;
    end
    Kvalues = KvaluesUnique;
    clear spIdxRow spIdxCol KvaluesUnique
    if buildMassMatrix
        Mvalues = MvaluesUnique;
        clear MvaluesUnique
    end
    spIdx = [spIdxRowUnique, spIdxColUnique];
    clear spIdxRowUnique spIdxColUnique
    
    [spIdx,~,I] = unique(spIdx,'rows','stable');
    Kvalues = accumarray(I,Kvalues);
    if buildMassMatrix
        Mvalues = accumarray(I,Mvalues);
    end
else
    spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
    spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
    Kvalues = reshape(Kvalues,numel(Kvalues),1);
    if buildMassMatrix
        Mvalues = reshape(Mvalues,numel(Mvalues),1);
    end
    spIdx = [spIdxRow, spIdxCol];
    clear spIdxRow spIdxCol
    [spIdx,~,I] = unique(spIdx,'rows','stable');
    Kvalues = accumarray(I,Kvalues);
    if buildMassMatrix
        Mvalues = accumarray(I,Mvalues);
    end
end

K = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Kvalues,noDofs,noDofs,numel(Kvalues));
clear Kvalues

if buildMassMatrix
    M = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Mvalues,noDofs,noDofs,numel(Mvalues));
else
    M = [];
end


