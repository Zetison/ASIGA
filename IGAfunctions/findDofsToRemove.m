function varCol = findDofsToRemove(varCol,discriminateIdenticalPoints)

if nargin < 2
    discriminateIdenticalPoints = false;
end
patches = varCol.patches;
noPatches = numel(patches);
noElemsPatch = zeros(noPatches,1);
noCtrlPtsPatch = zeros(noPatches,1);
d_p = numel(patches{1}.elRange);
noEl = zeros(noPatches,d_p);
for i = 1:noPatches
    noElemsPatch(i) = patches{i}.noElems;
    noCtrlPtsPatch(i) = patches{i}.noCtrlPts;
    for j = 1:d_p
        noEl(i,j) = size(patches{i}.elRange{j},1);
    end
end

noCtrlPts = sum(noCtrlPtsPatch);
noElems = sum(noElemsPatch);
pIndex = zeros(noElems,1);
isPML = false(noElems,1);
element = zeros(noElems,size(patches{1}.element,2));
controlPts = zeros(noCtrlPts,size(patches{1}.controlPts,2));
weights = zeros(noCtrlPts,1);
index = zeros(noElems,d_p);
elRange = cell(1,d_p);
for j = 1:d_p
    elRange{j} = zeros(sum(noEl(:,j)),2);
end
e = 1;
jC = 1;
maxDof = 0;
jEl = zeros(1,d_p);
knotVecs = cell(1,noPatches);
for i = 1:noPatches
    if isfield(varCol.nurbs{i},'isPML')
        isPML(e:e+noElemsPatch(i)-1) = varCol.nurbs{i}.isPML;
    end
    pIndex(e:e+noElemsPatch(i)-1) = i;
    element(e:e+noElemsPatch(i)-1,:) = maxDof + patches{i}.element;
    controlPts(jC:jC+noCtrlPtsPatch(i)-1,:) = patches{i}.controlPts;
    weights(jC:jC+noCtrlPtsPatch(i)-1,:) = patches{i}.weights;
    maxDof = maxDof + noCtrlPtsPatch(i);
    index(e:e+noElemsPatch(i)-1,:) = patches{i}.index + repmat(jEl,noElemsPatch(i),1);       
    for j = 1:d_p
        elRange{j}(jEl(j)+1:jEl(j)+noEl(i,j),:) = patches{i}.elRange{j};
    end
    e = e + noElemsPatch(i);
    jC = jC + noCtrlPtsPatch(i);
    jEl = jEl + noEl(i,:);

    knotVecs{i} = patches{i}.nurbs.knots;
end
d = varCol.dimension;
nodesMap = 1:noCtrlPts;
element2 = element;
if discriminateIdenticalPoints
    dofsToRemove = [];
    childrenNodes = [];
    gluedNodes = {};
else
    % Eps = 1e10*eps;
    Eps = 1e7*eps;
    % Eps = 1e5*eps;
    [~, gluedNodes] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(norm2(controlPts)), 'OutputAllIndices', true);
    repeatedNode = zeros(numel(gluedNodes),1);
    for i = 1:numel(gluedNodes)
        repeatedNode(i) = numel(gluedNodes{i}) - 1;
    end
    gluedNodes(repeatedNode == 0) = [];
    noChildrenNodes = sum(repeatedNode);
    childrenNodes = zeros(1,noChildrenNodes);

    counter = 1;
    if false % use slow method
        for i = 1:length(gluedNodes)
            parentIdx = gluedNodes{i}(1);
            for j = 2:length(gluedNodes{i})
                childrenIdx = gluedNodes{i}(j);
                indices = (element == childrenIdx);
                element(indices) = parentIdx;
                nodesMap(nodesMap == childrenIdx) = parentIdx;
                childrenNodes(counter) = childrenIdx;        
                counter = counter + 1;
            end
        end
    else % use fast method
        for i = 1:length(gluedNodes)
            childrenNodes_i = gluedNodes{i}(2:end);
            noChildrenNodes_i = numel(childrenNodes_i);
            nodesMap(childrenNodes_i) = gluedNodes{i}(1);
            childrenNodes(counter:counter+noChildrenNodes_i-1) = childrenNodes_i;
            counter = counter + noChildrenNodes_i;
        end
        element = nodesMap(element);
    end

    dofsToRemove = zeros(1,length(childrenNodes)*d);
    for i = 1:d
        dofsToRemove(i:d:end) = d*(childrenNodes-1)+i;
    end

    dofsToRemove = sort(unique(dofsToRemove));
end

varCol.childrenNodes = unique(childrenNodes);
varCol.noElemsPatch = noElemsPatch;
varCol.noCtrlPtsPatch = noCtrlPtsPatch;
varCol.noDofs = d*sum(noCtrlPtsPatch);
varCol.controlPts = controlPts;
varCol.weights = weights;
varCol.elRange = elRange;
varCol.index = index;
varCol.nodesMap = nodesMap;
varCol.element = element;
varCol.element2 = element2;
varCol.gluedNodes = gluedNodes;
varCol.degree = patches{1}.nurbs.degree;  % assume polynomial orders to be equal in all patches
varCol.knotVecs = knotVecs;
varCol.noCtrlPts = noCtrlPts;
varCol.noElems = noElems;
varCol.pIndex = pIndex;
varCol.noPatches = noPatches;
varCol.isPML = isPML;

if isfield(varCol,'geometry')
    varColBdry = meshBoundary(varCol,'homDirichlet');
    dirichletNodes = varColBdry.nodes;
    dofsToRemove = sort(unique([dofsToRemove,dirichletNodes]));
end
varCol.dofsToRemove = dofsToRemove;

