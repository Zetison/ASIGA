function varCol = findDofsToRemove(varCol)

patches = varCol.patches;
noPatches = numel(patches);
noElemsPatch = zeros(noPatches,1);
noCtrlPtsPatch = zeros(noPatches,1);
noParams = numel(patches{1}.elRange);
noEl = zeros(noPatches,noParams);
for i = 1:noPatches
    noElemsPatch(i) = patches{i}.noElems;
    noCtrlPtsPatch(i) = patches{i}.noCtrlPts;
    for j = 1:noParams
        noEl(i,j) = size(patches{i}.elRange{j},1);
    end
end


noCtrlPts = sum(noCtrlPtsPatch);
noElems = sum(noElemsPatch);
pIndex = zeros(noElems,1);
element = zeros(noElems,size(patches{1}.element,2));
controlPts = zeros(noCtrlPts,size(patches{1}.controlPts,2));
weights = zeros(noCtrlPts,1);
index = zeros(noElems,noParams);
elRange = cell(1,noParams);
for j = 1:noParams
    elRange{j} = zeros(sum(noEl(:,j)),2);
end
e = 1;
jC = 1;
maxDof = 0;
jEl = zeros(1,noParams);
knotVecs = cell(1,noPatches);
for i = 1:noPatches
    pIndex(e:e+noElemsPatch(i)-1) = i;
    element(e:e+noElemsPatch(i)-1,:) = maxDof + patches{i}.element;
    controlPts(jC:jC+noCtrlPtsPatch(i)-1,:) = patches{i}.controlPts;
    weights(jC:jC+noCtrlPtsPatch(i)-1,:) = patches{i}.weights;
    maxDof = maxDof + noCtrlPtsPatch(i);
    index(e:e+noElemsPatch(i)-1,:) = patches{i}.index + repmat(jEl,noElemsPatch(i),1);       
    for j = 1:noParams
        elRange{j}(jEl(j)+1:jEl(j)+noEl(i,j),:) = patches{i}.elRange{j};
    end
    e = e + noElemsPatch(i);
    jC = jC + noCtrlPtsPatch(i);
    jEl = jEl + noEl(i,:);

    knotVecs{i} = patches{i}.nurbs.knots;
end
d = varCol.dimension;
% Eps = 1e7*eps;
Eps = 1e5*eps;
% [~, I, IC] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(norm2(controlPts)));
% [~, I, IC] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(max(abs(controlPts))));
[~, gluedNodes] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(norm2(controlPts)), 'OutputAllIndices', true);
repeatedNode = zeros(numel(gluedNodes),1);
for i = 1:numel(gluedNodes)
    repeatedNode(i) = numel(gluedNodes{i}) - 1;
end
gluedNodes(repeatedNode == 0) = [];
noChildrenNodes = sum(repeatedNode);
% nI = setdiff(1:noCtrlPts,I);
% Iunique = I(IC);
% Im = unique(Iunique(nI));
% gluedNodes = cell(length(Im),1);
% 
% parfor i = 1:length(Im)
%     temp = repmat(controlPts(Im(i),:),noCtrlPts,1);
%     temp = controlPts-temp;
%     gluedNodes{i} = find(norm2(temp)./norm2(controlPts) < Eps);
% end
% childrenNodes = zeros(size(nI));
childrenNodes = zeros(1,noChildrenNodes);
element2 = element;
counter = 1;
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        childrenIdx = gluedNodes{i}(j);
        indices = (element == childrenIdx);
        element(indices) = parentIdx;
        childrenNodes(counter) = childrenIdx;        
        counter = counter + 1;
    end
end

dofsToRemove = zeros(1,length(childrenNodes)*d);
for i = 1:d
    dofsToRemove(i:d:end) = d*(childrenNodes-1)+i;
end

dofsToRemove = sort(unique(dofsToRemove));

varCol.noElemsPatch = noElemsPatch;
varCol.noCtrlPtsPatch = noCtrlPtsPatch;
varCol.noDofs = d*sum(noCtrlPtsPatch);
varCol.controlPts = controlPts;
varCol.weights = weights;
varCol.elRange = elRange;
varCol.index = index;
varCol.element = element;
varCol.element2 = element2;
varCol.dofsToRemove = dofsToRemove;
varCol.gluedNodes = gluedNodes;
varCol.degree = patches{1}.nurbs.degree;  % assume polynomial orders to be equal in all patches
varCol.knotVecs = knotVecs;
varCol.noCtrlPts = noCtrlPts;
varCol.noElems = noElems;
varCol.pIndex = pIndex;
varCol.noPatches = noPatches;

