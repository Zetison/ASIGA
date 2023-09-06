function relError = calcSurfErrorVec(task, i_domain, i_o)

LpOrder = task.err.LpOrder;
degree = task.varCol{i_domain}.degree(1:2);
weights = task.varCol{i_domain}.weights;
controlPts = task.varCol{i_domain}.controlPts;

if nargin < 2
    U = task.varCol{i_domain}.U;
else
    U = task.varCol{i_domain}.U(:,i_o);
end

varColBdry = meshBoundary(task.varCol{i_domain},'Gamma');

zeta0Nodes = varColBdry.nodes;
noElems = varColBdry.noElems;
element = varColBdry.element;
element2 = varColBdry.element2;
index = varColBdry.index;
pIndex = varColBdry.pIndex;
knotVecs = varColBdry.knotVecs;
elRange = varColBdry.elRange;

extraGP = task.misc.extraGP(1:2);
[Q, W] = gaussTensorQuad(degree(1:2)+3+extraGP);

p_h = complex(zeros(size(W,1),noElems,size(U,2)));
fact = zeros(size(W,1),noElems);
points = zeros(size(W,1), noElems,3);

for e = 1:noElems
% parfor e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(2,2);
    for i = 1:2
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = zeta0Nodes(element(e,:));
    
    pts = controlPts(sctr,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:); % New
    
    U_sctr = U(sctr,:);

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    
    J_1 = getJacobian(R,pts,2);
    
    p_h(:,e,:) = R{1}*U_sctr;
    fact(:,e) = J_1 * J_2 .* W;
    points(:,e,:) = R{1}*pts;
end
p_h = reshape(p_h, size(p_h,1)*size(p_h,2),size(U,2));
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
analyticFunctions = task.analyticFunctions({points});
p = analyticFunctions{1}.p;

if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*fact)).^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*fact)).^(1/LpOrder);
end

relError = 100*Error./normalization;
