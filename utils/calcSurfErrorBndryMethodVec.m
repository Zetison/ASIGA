function relError = calcSurfErrorBndryMethodVec(varCol, LpOrder)

i = 1;
degree = varCol{i}.degree; % assume degree is equal in all patches

index = varCol{i}.index;
noElems = varCol{i}.noElems;
elRange = varCol{i}.elRange;
element = varCol{i}.element;
element2 = varCol{i}.element2;
weights = varCol{i}.weights;
controlPts = varCol{i}.controlPts;
knotVecs = varCol{i}.knotVecs;
pIndex = varCol{i}.pIndex;
d_p = varCol{i}.patches{1}.nurbs.d_p;
U = varCol{i}.U;
k = varCol{i}.k;

if strcmp(varCol{1}.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol{1}.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end

extraGP = varCol{1}.extraGP;
[Q, W] = gaussTensorQuad(degree+3+extraGP);

p_h = complex(zeros(size(W,1), noElems, size(U,2)));
fact = zeros(size(W,1), noElems);
points = zeros(size(W,1), noElems,3);
isMFS = strcmp(varCol{1}.method,'MFS');

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
    y_s = varCol{i}.y_s;
    r = zeros(size(points));
    for j = 1:n_cp
        r(:,j) = norm2(elementAddition(-y_s(j,:), points));
    end
    p_h = Phi_k(r)*U;
else
    p_h = reshape(p_h, size(p_h,1)*size(p_h,2),size(U,2));
end

if varCol{i}.solveForPtot && varCol{i}.exteriorProblem
    analytic = @(x) varCol{i}.p_(x) + varCol{i}.p_inc(x);
else
    analytic = varCol{i}.p_;
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


