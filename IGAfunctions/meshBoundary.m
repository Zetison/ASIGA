function varColBdry = meshBoundary(varCol,name)
if ~isa(name,'char')
    error('This parameter should be a string')
end
gluedNodes = varCol.gluedNodes;

degree = varCol.degree(1:2);
patches = varCol.patches;

p_xi = degree(1); % assume p_xi is equal in all patches
p_eta = degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);


noParams = 2;
noElemsPatch = varCol.noElemsPatch;
set = varCol.geometry.topologysets.set;
noElems = varCol.noElems;
pIndex = zeros(noElems,1);
element = zeros(noElems,n_en);
index = zeros(noElems,noParams);
eBdry = 1;
maxDof = 0;
jEl = zeros(1,2);
elemMap = zeros(noElems,1);
nodes = zeros(1,varCol.noDofs);
counter = 1;
[idx,setFound] = findSet(set,name);
if ~setFound
    varColBdry.nodes = [];
    return
end
noItems = numel(set{idx}.item);
sub_nurbs = cell(1,noItems);

for i_free = 1:noItems
    patch = set{idx}.item{i_free}.Attributes.patch;
    midx = set{idx}.item{i_free}.Text;
    nurbs = patches{patch}.nurbs;
    bdryNodes = extractBdryNodes(nurbs,midx) + sum(varCol.noCtrlPtsPatch(1:patch-1));
    nodes(counter:counter+numel(bdryNodes)-1) = bdryNodes;
    at = zeros(2,3,'logical');
    at(midx) = true;
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = subNURBS({nurbs},'at',at.');
    sub_nurbs(i_free) = varCol_dummy.nurbs;
    varCol_dummy = generateIGAmesh(convertNURBS(varCol_dummy));
    noElemsXiEta = varCol_dummy.patches{1}.noElems;
    index(eBdry:eBdry+noElemsXiEta-1,:) = varCol_dummy.patches{1}.index + repmat(jEl,noElemsXiEta,1);    
    pIndex(eBdry:eBdry+noElemsXiEta-1) = patch;
    element(eBdry:eBdry+noElemsXiEta-1,:) = maxDof + varCol_dummy.patches{1}.element;
    noEl1 = size(varCol_dummy.patches{1}.elRange{1},1);
    noEl2 = size(varCol_dummy.patches{1}.elRange{2},1);
    noEl3 = size(patches{patch}.elRange{ceil(midx/2)},1);
    jEl = jEl + [noEl1,noEl2];
    maxDof = maxDof + numel(bdryNodes);
    
    temp = repmat(eBdry:eBdry+noElemsXiEta-1,1,noEl3);
    e = sum(noElemsPatch(1:patch-1)) + 1;    
    elemMap(e:e+sum(noElemsPatch(patch))-1) = temp(:);
    
    counter = counter+numel(bdryNodes);
    eBdry = eBdry + noElemsXiEta;
end



nodes(counter:end) = [];
pIndex(eBdry:end) = [];
element(eBdry:end,:) = [];
index(eBdry:end,:) = [];
noElemsBdry = numel(pIndex);
noSurfDofs = numel(nodes);

element2 = element;
% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    indices = (nodes(element(:)) == gluedNodes{i}(1));
    if any(indices)
        parentIdx = element(indices);
        for j = 2:length(gluedNodes{i})
            indices = (nodes(element(:)) == gluedNodes{i}(j));
            element(indices) = parentIdx(1);
        end
    end
end
varColBdry.nodes = nodes;
varColBdry.noElems = noElemsBdry;
varColBdry.element = element;
varColBdry.element2 = element2;
varColBdry.index = index;
varColBdry.pIndex = pIndex;
varColBdry.n_en = n_en;
varColBdry.noSurfDofs = noSurfDofs;
varColBdry.elemMap = elemMap;
varColBdry.sub_nurbs = sub_nurbs;

function bdryNodes = extractBdryNodes(nurbs,midx)
n_xi = nurbs.number(1);
n_eta = nurbs.number(2);
n_zeta = nurbs.number(3);
noControlPts = n_xi*n_eta*n_zeta;
nodes = reshape(1:noControlPts,n_xi,n_eta,n_zeta);
switch midx
    case 1
        bdryNodes = nodes(1,:,:);
    case 2
        bdryNodes = nodes(end,:,:);
    case 3
        bdryNodes = nodes(:,1,:);
    case 4
        bdryNodes = nodes(:,end,:);
    case 5
        bdryNodes = nodes(:,:,1);
    case 6
        bdryNodes = nodes(:,:,end);
end
bdryNodes = bdryNodes(:);






