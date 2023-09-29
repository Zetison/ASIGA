function I_tot = integrateFunc(task,fun,fun_ref)

weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;
index = task.varCol{1}.index;
elRange = task.varCol{1}.elRange;
noElems = task.varCol{1}.noElems;
degree = task.varCol{1}.degree;
knotVecs = task.varCol{1}.knotVecs;
pIndex = task.varCol{1}.pIndex;
element = task.varCol{1}.element;
element2 = task.varCol{1}.element2;
set = task.varCol{1}.geometry.topologysets.set;
setIdx = findSet(set,'Gamma_b');

item = set{setIdx}.item;
noPMLpatches = numel(item);
Gamma_b_patches = zeros(noPMLpatches,1);
for i = 1:noPMLpatches
    Gamma_b_patches(i) = item{i}.Attributes.patch;
end
surfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRange{3}(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(1) == 0 && ismember(pIndex(e),Gamma_b_patches)
        surfaceElements = [surfaceElements e];
    end
end

% varColBdry = meshBoundary(task.varCol{1},'Gamma_b');
    
% zeta0Nodes = varColBdry.nodes;
% element = varColBdry.element;
% element2 = varColBdry.element2;
% index = varColBdry.index;
% pIndex = varColBdry.pIndex;
% elRange = varColBdry.elRange;

extraGP = task.misc.extraGP;
noGp = degree+3+extraGP(1:numel(degree));
[Q, W] = gaussTensorQuad(noGp(1:2)); 

Error = 0;
normalization = 0;
% for i = 1:length(surfaceElements) %
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

    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;

    xi = [parent2ParametricSpace(Xi_e, Q), ones(size(Q,1),1)];
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1,crossProd] = getJacobian(R,pts,2);

    n = crossProd./repmat(J_1,1,3);
    if any(abs(norm2(n) - 1) > 1e4*eps)
        keyboard
    end
    y = R{1}*pts;  
    dXdzeta = R{4}*pts;
    Error = Error + (norm2(fun(y, n, dXdzeta) - fun_ref(y, n, dXdzeta)).^2 .*J_1 * J_2).' * W;
    normalization = normalization + (norm2(fun_ref(y, n, dXdzeta)).^2 .*J_1 * J_2).' * W;
end

I_tot = 100*sqrt(Error./normalization);


