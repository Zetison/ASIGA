function task = applyNeumannCondition(task,i_domain)
noDofs = task.varCol{i_domain}.noDofs;
weights = task.varCol{i_domain}.weights;
controlPts = task.varCol{i_domain}.controlPts;
degree = task.varCol{i_domain}.degree(1:2);
elRange = task.varCol{i_domain}.elRange;
knotVecs = task.varCol{i_domain}.knotVecs;
d_f = task.varCol{i_domain}.fieldDimension;
media = task.varCol{i_domain}.media;
useROM = task.rom.useROM;
noRHSs = task.noRHSs;

switch media
    case 'fluid' 
        if useROM
            dp_inc = task.dp_inc_ROM_;  
        else
            dp_inc = task.dp_inc_;  
        end
        p_inc = NaN; % for for-loop transparency
    case 'solid'
        if useROM
            p_inc = task.dp_inc_;  
        else
            p_inc = task.p_inc_;  
        end
        dp_inc = NaN; % for for-loop transparency
end

if task.varCol{i_domain}.boundaryMethod
    noElems = task.varCol{i_domain}.noElems;
    element = task.varCol{i_domain}.element;
    element2 = task.varCol{i_domain}.element2;
    index = task.varCol{i_domain}.index;
    pIndex = task.varCol{i_domain}.pIndex;
    n_en = prod(degree+1);
    zeta0Nodes = 1:noDofs;
else
    [zeta0Nodes, noElems, element, element2, index, pIndex, n_en] = meshBoundary(task.varCol{i_domain},'Neumann');
end
Fvalues = zeros(d_f*n_en,noElems,noRHSs);
indices = zeros(d_f*n_en,noElems);

extraGP = task.misc.extraGP;
[Q, W] = gaussTensorQuad(degree+1+extraGP(1:2));

parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(2,2);
    for i = 1:2
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = zeta0Nodes(element(e,:));
    sctr_f_e = zeros(d_f*n_en,1);
    for i = 1:d_f
        sctr_f_e(i:d_f:d_f*n_en) = d_f*(sctr-1)+i;
    end
    
    pts = controlPts(sctr,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:); % New
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,2);
    fact = J_1 * J_2 .* W;
    X = R{1}*pts;
    n = crossProd./repmat(J_1,1,3);
    switch media
        case 'fluid' 
            Fvalues(:,e,:) = R{1}.'*(dp_inc(X,n).*fact);
        case 'solid'
            Fvalues(:,e,:) = -kron2(n,R{1})*(p_inc(X).*fact);
    end
    
    indices(:,e) = sctr_f_e.'; 
end

F = zeros(noDofs,noRHSs);
for i = 1:noRHSs
    F(:,i) = vectorAssembly(Fvalues(:,:,i),indices,noDofs);
end
task.varCol{i_domain}.FF = F;



