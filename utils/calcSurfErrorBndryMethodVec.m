function relError = calcSurfErrorBndryMethodVec(task)

i = 1;
LpOrder = task.err.LpOrder;
degree = task.varCol{i}.degree; % assume degree is equal in all patches

index = task.varCol{i}.index;
noElems = task.varCol{i}.noElems;
elRange = task.varCol{i}.elRange;
element = task.varCol{i}.element;
element2 = task.varCol{i}.element2;
weights = task.varCol{i}.weights;
controlPts = task.varCol{i}.controlPts;
knotVecs = task.varCol{i}.knotVecs;
pIndex = task.varCol{i}.pIndex;
d_p = task.varCol{i}.patches{1}.nurbs.d_p;
U = task.varCol{i}.U;
k = task.misc.omega/task.varCol{i}.c_f;

if strcmp(task.misc.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = task.varCol{1}.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end

extraGP = task.misc.extraGP;
[Q, W] = gaussTensorQuad(degree+3+extraGP(1:2));

p_h = complex(zeros(size(W,1), noElems, size(U,2)));
fact = zeros(size(W,1), noElems);
points = zeros(size(W,1), noElems,3);
isMFS = strcmp(task.misc.method,'MFS');

% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(d_p,2);
    for ii = 1:d_p
        Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;


    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    J_1 = getJacobian(R,pts,d_p);
        
    y = R{1}*pts;
    if useEnrichedBfuns
        R{1} = R{1}*repmat(exp(1i*k*(y*d_vec)),1,size(R{1},2));
    end
    if ~isMFS
        U_sctr = U(sctr,:);
        p_h(:,e,:) = R{1}*U_sctr;
    end
    fact(:,e) = J_1 * J_2 .* W;
    points(:,e,:) = y;
end
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
if isMFS
    Phi_k = @(r) exp(1i*k*r)./(4*pi*r);
    n_cp = size(U,1);
    y_s = task.varCol{i}.y_s;
    r = zeros(size(points));
    for j = 1:n_cp
        r(:,j) = norm2(elementAddition(-y_s(j,:), points));
    end
    p_h = Phi_k(r)*U;
else
    p_h = reshape(p_h, size(p_h,1)*size(p_h,2),size(U,2));
end

if task.misc.solveForPtot && task.misc.exteriorProblem
    analytic = @(x) task.varCol{i}.p_(x) + task.p_inc_(x);
else
    analytic = task.varCol{i}.p_;
end
p = analytic(points);
if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*fact)).^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*fact)).^(1/LpOrder);
end
relError = 100*Error./normalization;


