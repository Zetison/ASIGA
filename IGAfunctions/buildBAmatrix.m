function [M, F, varCol] = buildBAmatrix(varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.


%% Extract all needed data from options and varCol
degree = varCol.degree; % assume degree is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
noDofs = varCol.noCtrlPts;
extraGP = varCol.extraGP;
d_p = varCol.patches{1}.nurbs.d_p;

if varCol.solveForPtot && varCol.exteriorProblem
    analytic = @(x) varCol.p(x) + varCol.p_inc(x);
else
    analytic = varCol.p;
end

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
k = varCol.k;
totNoQP = 0;
[Q, W] = gaussTensorQuad(degree+1+extraGP);
n_en = prod(degree+1);

%% Preallocation and initiallizations
v_values   = zeros(size(W,1),noElems,3); 
n_values   = zeros(size(W,1),noElems,3); 
J_1_values   = zeros(size(W,1),noElems,1); 

progressBars = varCol.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    ppm = ParforProgMon('Building BA matrix (1/2): ', noElems, nProgressStepSize);
else
    ppm = NaN;
end


%% Find nodes at which to evaluate the exact solution
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

    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1,crossProd] = getJacobian(R,pts,d_p);
    n_values(:,e,:) = crossProd./J_1;
    v_values(:,e,:) = R{1}*pts;
    J_1_values(:,e) = J_1;
end
v_values = reshape(v_values, size(W,1)*noElems, 3);
n_values = reshape(n_values, size(W,1)*noElems, 3);
if nargin(analytic) == 2
    analytic_values = analytic(v_values,n_values);
else
    analytic_values = analytic(v_values);
end

no_funcs = size(analytic_values,2);
analytic_values = permute(reshape(analytic_values, size(W,1), noElems, no_funcs),[1,3,2]);


%% Build global matrices
sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
Fvalues   = zeros(n_en,noElems,no_funcs); 
F_indices = zeros(n_en,noElems); 
if progressBars
    ppm = ParforProgMon('Building BA matrix (2/2): ', noElems, nProgressStepSize);
else
    ppm = NaN;
end

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
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);

    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    
    J_1 = J_1_values(:,e);

    y = R{1}*pts;
    if useEnrichedBfuns
        temp = exp(1i*k*(y*d_vec));
        R{1} = R{1}.*temp(:,ones(1,noGp));
    end
    m_e = zeros(n_en);
    for i = 1:numel(W)
        m_e = m_e + R{1}(i,:)'*R{1}(i,:) * abs(J_1(i)) * J_2 * W(i);  
    end

    spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
    spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
    Mvalues(:,e) = reshape(m_e, sizeMe, 1);

    F_indices(:,e) = sctr';
    Fvalues(:,e,:) = R{1}.'*(analytic_values(:,:,e).*repmat(abs(J_1)*J_2.*W,1,no_funcs));
    totNoQP = totNoQP + size(Q,1);
end
        

%% Collect data into global matrices (and load vector)
F = vectorAssembly(Fvalues,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);

varCol.totNoQP = totNoQP;

