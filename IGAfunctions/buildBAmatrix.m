function task = buildBAmatrix(task,i_varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.


%% Extract all needed data from options and varCol
degree = task.varCol{i_varCol}.degree; % assume degree is equal in all patches

index = task.varCol{i_varCol}.index;
noElems = task.varCol{i_varCol}.noElems;
elRange = task.varCol{i_varCol}.elRange;
element = task.varCol{i_varCol}.element;
element2 = task.varCol{i_varCol}.element2;
weights = task.varCol{i_varCol}.weights;
controlPts = task.varCol{i_varCol}.controlPts;
knotVecs = task.varCol{i_varCol}.knotVecs;
pIndex = task.varCol{i_varCol}.pIndex;
noDofs = task.varCol{i_varCol}.noDofs;
extraGP = task.misc.extraGP;
d_p = task.varCol{i_varCol}.patches{1}.nurbs.d_p;

if task.misc.solveForPtot && task.misc.exteriorProblem
    analytic = @(x) task.varCol{i_varCol}.p_(x) + task.p_inc_(x);
else
    switch task.varCol{i_varCol}.media
        case 'fluid'
            analytic = task.varCol{i_varCol}.p_;
        case 'solid'
            analytic = @(X) reshape([task.varCol{i_varCol}.u_x_(X).'; task.varCol{i_varCol}.u_y_(X).'; task.varCol{i_varCol}.u_z_(X).'],3*size(X{i_varCol},1),1);
    end
end

if strcmp(task.misc.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = task.varCol{i_varCol}.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
if strcmp(task.varCol{i_varCol}.media,'fluid')
    k = task.misc.omega/task.varCol{i_varCol}.c_f;
else
    k = NaN;
end
[Q, W] = gaussTensorQuad(degree+1+extraGP(1:d_p));
n_en = prod(degree+1);

%% Preallocation and initiallizations
v_values   = zeros(size(W,1),noElems,3); 
n_values   = zeros(size(W,1),noElems,3); 
J_1_values   = zeros(size(W,1),noElems,1); 

progressBars = task.misc.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    try
        ppm = ParforProgMon('Building BA matrix (1/2): ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
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
X = cell(1,i_varCol);
for i = 1:i_varCol
    if i == i_varCol
        X{i} = v_values;
    else
        X{i} = zeros(0,3);
    end
end
if nargin(analytic) == 2
    analytic_values = analytic(X,n_values);
else
    analytic_values = analytic(X);
end

d_f = task.varCol{i_varCol}.dimension;
analytic_values = reshape(analytic_values, d_f, size(W,1), noElems);


%% Build RHS
Fvalues   = zeros(n_en*d_f,noElems); 
F_indices = zeros(n_en*d_f,noElems); 
if progressBars
    try
        ppm = ParforProgMon('Building BA matrix (2/2): ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
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
    
    J_1 = J_1_values(:,e);

    y = R{1}*pts;
    if useEnrichedBfuns
        temp = exp(1i*k*(y*d_vec));
        R{1} = R{1}.*temp(:,ones(1,noGp));
    end

    F_indices(:,e) = sctr_k_e.';
    Fvalues(:,e) = reshape(analytic_values(:,:,e)*(abs(J_1)*J_2.*W.*R{1}),n_en*d_f,1);
end
        

%% Collect data into global matrices (and load vector)
task.varCol{i_varCol}.FF = vectorAssembly(Fvalues,F_indices,noDofs);

