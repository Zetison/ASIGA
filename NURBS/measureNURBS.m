function [A, Ae] = measureNURBS(nurbs)
varCol.dimension = 1;
varCol.nurbs = nurbs;
varCol = convertNURBS(varCol);
varCol = generateIGAmesh(varCol);
varCol = findDofsToRemove(varCol);

%% Extract all needed data from varCol
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

%% Preallocation and initiallizations
[Q,W] = gaussTensorQuad(degree+20); 
d_p = varCol.patches{1}.nurbs.d_p;
Ae = zeros(noElems,1);
% for e = 1:noElems
parfor e = 1:noElems
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
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    J_1 = getJacobian(R,pts,d_p);
    Ae(e) = J_1.' * J_2 * W;
end
A = sum(Ae);
