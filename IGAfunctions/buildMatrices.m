function varCol = buildMatrices(varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles elasticity, laplace- and poisson equation, 
% and dynamic versions of these.

%% Extract all needed data from options and varCol
degree = varCol.degree;
knotVecs = varCol.knotVecs;
index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
pIndex = varCol.pIndex;
noDofs = varCol.noDofs;
noCtrlPts = varCol.noCtrlPts;

d_f = varCol.fieldDimension;
d_p = varCol.patches{1}.nurbs.d_p;
buildMassMatrix = varCol.buildMassMatrix;
applyBodyLoading = varCol.applyBodyLoading;
operator = varCol.operator;

if strcmp(operator,'linearElasticity')
    C = varCol.C;
else
    C = NaN; % Will not be used
end

%% Preallocation and initiallizations
n_en = prod(degree+1);
sizeKe = (d_f*n_en)^2;
sizeMe = n_en^2;
spIdxRow = zeros(sizeKe,noElems,'uint32');
spIdxCol = zeros(sizeKe,noElems,'uint32');
spIdxRowM = zeros(sizeMe,noElems,'uint32');
spIdxColM = zeros(sizeMe,noElems,'uint32');
Kvalues = zeros(sizeKe,noElems); 

if buildMassMatrix
    Mvalues = zeros(sizeMe,noElems); 
end
if applyBodyLoading
    F_indices = zeros(d_f*n_en,noElems); 
    Fvalues = zeros(d_f*n_en,noElems); 
    force = varCol.force;
else
    force = NaN;
end

extraGP = varCol.extraGP;
[Q, W] = gaussTensorQuad(degree+1+extraGP);

progressBars = varCol.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    ppm = ParforProgMon('Building mass/stiffness matrix: ', noElems, nProgressStepSize);
else
    ppm = NaN;
end

%% Build global matrices
% for e = 1:noElems
parfor e = 1:noElems
	if progressBars && mod(e,nProgressStepSize) == 0
        ppm.increment();
	end
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(d_p,2);
    for i = 1:d_p
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;
    
    sctr_k_e = zeros(1,d_f*n_en);
    for i = 1:d_f
        sctr_k_e(i:d_f:end) = d_f*(sctr-1)+i;
    end

    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    J_1 = getJacobian(R,pts,d_p);
    fact = J_1 * J_2 .* W;

    dXdxi = R{2}*pts;
    dXdeta = R{3}*pts;
    dXdzeta = R{4}*pts;
                
    a11 = dXdxi(:,1);
    a21 = dXdxi(:,2);
    a31 = dXdxi(:,3);
    a12 = dXdeta(:,1);
    a22 = dXdeta(:,2);
    a32 = dXdeta(:,3);
    a13 = dXdzeta(:,1);
    a23 = dXdzeta(:,2);
    a33 = dXdzeta(:,3);
    Jinv1 = [(a22.*a33-a23.*a32)./J_1, (a23.*a31-a21.*a33)./J_1, (a21.*a32-a22.*a31)./J_1];
    Jinv2 = [(a13.*a32-a12.*a33)./J_1, (a11.*a33-a13.*a31)./J_1, (a12.*a31-a11.*a32)./J_1];
    Jinv3 = [(a12.*a23-a13.*a22)./J_1, (a13.*a21-a11.*a23)./J_1, (a11.*a22-a12.*a21)./J_1];
    dRdX = cell(d_f,1);
    dRdX{1} = repmat(Jinv1(:,1),1,n_en).*R{2} + repmat(Jinv1(:,2),1,n_en).*R{3} + repmat(Jinv1(:,3),1,n_en).*R{4};
    dRdX{2} = repmat(Jinv2(:,1),1,n_en).*R{2} + repmat(Jinv2(:,2),1,n_en).*R{3} + repmat(Jinv2(:,3),1,n_en).*R{4};
    dRdX{3} = repmat(Jinv3(:,1),1,n_en).*R{2} + repmat(Jinv3(:,2),1,n_en).*R{3} + repmat(Jinv3(:,3),1,n_en).*R{4};
    
    Kvalues(:,e) = stiffnessElementMatrix(dRdX,fact,d_f,n_en,operator,C);
    if buildMassMatrix
        Mvalues(:,e) = kron2(R{1},R{1}) * fact;
    end
    if applyBodyLoading
        v = R{1}*pts;
        f_gp = force(v);
        F_indices(:,e) = sctr_k_e;
        Fvalues(:,e) = kron2(f_gp,R{1}) * fact;
    end
    spIdxRow(:,e) = kron(ones(1,d_f*n_en),sctr_k_e);
    spIdxCol(:,e) = kron(sctr_k_e,ones(1,d_f*n_en));
    spIdxRowM(:,e) = kron(ones(1,n_en),sctr);
    spIdxColM(:,e) = kron(sctr,ones(1,n_en));
end

%% Collect data into global matrices (and load vector)
varCol.A_K = sparse(double(spIdxRow),double(spIdxCol),Kvalues,noDofs,noDofs,numel(Kvalues));

if applyBodyLoading
    varCol.F = vectorAssembly(Fvalues,F_indices,noDofs);
end

if buildMassMatrix
    varCol.A_M = kron(sparse(double(spIdxRowM),double(spIdxColM),Mvalues,noCtrlPts,noCtrlPts,numel(Mvalues)),eye(d_f));
end

