function varCol = applyNeumannCondition(varCol,outerBoundary)
noDofs = varCol.noDofs;
weights = varCol.weights;
controlPts = varCol.controlPts;
degree = varCol.degree(1:2);
elRange = varCol.elRange;
d_p = varCol.patches{1}.nurbs.d_p;
knotVecs = varCol.knotVecs;
d_f = varCol.fieldDimension;
media = varCol.media;
useROM = varCol.useROM;
noRHSs = varCol.noRHSs;
if useROM
    dp_inc_ROM_ = varCol.dp_inc_ROM_;
    p_inc_ROM_ = varCol.p_inc_ROM_;
end
switch media
    case 'fluid' 
        if useROM
            dp_inc = @(X,n) dp_inc_ROM_(X,n);  
        else
            dp_inc = @(X,n) varCol.dp_inc_(X,n);  
        end
        p_inc = NaN; % for for-loop transparency
    case 'solid'
        if useROM
            p_inc = @(X) p_inc_ROM_(X);  
        else
            p_inc = @(X) varCol.p_inc_(X);  
        end
        dp_inc = NaN; % for for-loop transparency
end

if varCol.boundaryMethod
    noElems = varCol.noElems;
    element = varCol.element;
    element2 = varCol.element2;
    index = varCol.index;
    pIndex = varCol.pIndex;
    n_en = prod(degree+1);
    noSurfDofs = noDofs;
    noDofs = 0;
    zeta0Nodes = 1:noSurfDofs;
else
    [zeta0Nodes, noElems, element, element2, index, pIndex, n_en] = meshBoundary(varCol,outerBoundary);
end
Fvalues = zeros(d_f*n_en,noElems,noRHSs);
indices = zeros(d_f*n_en,noElems);

[Q, W] = gaussTensorQuad(degree+1);

parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(d_p-1,2);
    for i = 1:d_p-1
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = zeta0Nodes(element(e,:));
    sctr_f_e = zeros(d_f*n_en,1);
    for i = 1:d_f
        sctr_f_e(i:d_f:d_f*n_en) = d_f*(sctr-1)+i;
    end
    
    pts = controlPts(sctr,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:); % New
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^(d_p-1);
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,d_p-1);
    fact = J_1 * J_2 .* W;
    X = R{1}*pts;
    n = crossProd./repmat(J_1,1,3);
    switch media
        case 'fluid' 
            Fvalues(:,e,:) = R{1}'*(dp_inc(X,n).*fact);
        case 'solid'
            Fvalues(:,e,:) = -kron2(n,R{1})*(p_inc(X).*fact);
    end
    
    indices(:,e) = sctr_f_e.'; 
end

F = zeros(noDofs,noRHSs);
for i = 1:noRHSs
    F(:,i) = vectorAssembly(Fvalues(:,:,i),indices,noDofs);
end
varCol.FF = F;



