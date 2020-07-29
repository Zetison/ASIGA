function [K, M, v_arr] = buildGlobalMatrices_forBAanalysis(varCol, newOptions)
error('Depricated. Use buildBAmatrix instead')
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

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
l = varCol.nurbs.number(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;

noCtrlPts = varCol.noCtrlPts;
noDofs = varCol.noDofs;

if strcmp(options.operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = (p+1)*(q+1)*(r+1);

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

[W3D,Q3D] = gaussianQuadNURBS(p+1+floor(16*p/n),q+1+floor(8*q/m),r+1+floor(4*r/l)); 
% [W3D,Q3D] = gaussianQuadNURBS(p+1,q+1,r+1); 

v_arr = zeros(size(W3D,1), noElems, 3);

%% Build global matrices
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
    v_arr_temp = zeros(size(W3D,1), 3);
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
        
        switch options.operator
            case 'linearElasticity'
                B = strainDispMatrix3d(n_en,dRdX);
                k_e = k_e + B' * C * B * abs(J_1) * J_2 * wt; 
                if options.buildMassMatrix
                    m_e = m_e + blkdiag(R_fun'*R_fun, R_fun'*R_fun, R_fun'*R_fun) * abs(J_1) * J_2 * wt;      
                end
            case 'laplace'
                k_e = k_e + dRdX'*dRdX* abs(J_1) * J_2 * wt;  
                if options.buildMassMatrix
                    m_e = m_e + R_fun'*R_fun * abs(J_1) * J_2 * wt;  
                end
        end
        
        v_arr_temp(gp,:) = R_fun*pts;
        if options.applyBodyLoading
            v = R_fun*pts;
            f_gp = varCol.f(v(1),v(2),v(3));
            f_e = f_e + [f_gp(1)*R_fun'; f_gp(2)*R_fun'; f_gp(3)*R_fun'] * abs(J_1) * J_2 * wt;
        end
    end
    v_arr(:,e,:) = v_arr_temp;
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
v_arr = reshape(v_arr, size(v_arr,1)*size(v_arr,2),3);
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
