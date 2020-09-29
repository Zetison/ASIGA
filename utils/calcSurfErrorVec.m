function relError = calcSurfErrorVec(varCol, U, LpOrder)

degree = varCol.degree;

index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;


% Find elements on the inner surface for evaluation of
% backscattered pressure in far field
surfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRange{3}(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(1) == 0
        surfaceElements = [surfaceElements e];
    end
end

% if noElems < 600
%     noQuadPts = 16;
% else
%     noQuadPts = 10;
% end

extraGP = varCol.extraGP;
[Q, W] = gaussTensorQuad(degree(1:2)+3+extraGP);

p_h = zeros(size(W,1),length(surfaceElements));
fact = zeros(size(W,1),length(surfaceElements));
points = zeros(size(W,1), length(surfaceElements),3);

% for i = 1:length(surfaceElements)
parfor i = 1:length(surfaceElements)
    e = surfaceElements(i);
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(2,2);
    for ii = 1:2
        Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);
    U_sctr = U(sctr,:);

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;

    xi = [parent2ParametricSpace(Xi_e, Q), zeros(size(Q,1),1)];
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    J_1 = getJacobian(R,pts,2);
    
    p_h(:,i) = R{1}*U_sctr;
    fact(:,i) = J_1 * J_2 .* W;
    points(:,i,:) = R{1}*pts;
end
p_h = reshape(p_h, size(p_h,1)*size(p_h,2),1);
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
analyticFunctions = varCol.analyticFunctions({points});
p = analyticFunctions{1}.p;

if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*fact))^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*fact))^(1/LpOrder);
end

relError = 100*Error/normalization;


