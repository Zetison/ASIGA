function varColBdry = meshBoundary(varCol,name)
if ~isa(name,'char')
    error('This parameter should be a string')
end

degree = varCol.degree(1:2);
patches = varCol.patches;

p_xi = degree(1); % assume p_xi is equal in all patches
p_eta = degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);

d_p = 2;
noElemsPatch = varCol.noElemsPatch;
set = varCol.geometry.topologysets.set;
connection = varCol.geometry.topology.connection;
noElems = varCol.noElems;
pIndex = zeros(noElems,1);
element = zeros(noElems,n_en);
index = zeros(noElems,d_p);
eBdry = 1;
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
elRange = cell(1,d_p);
knotVecs = cell(1,noItems);
for i = 1:d_p
    elRange{i} = zeros(noElems,d_p);
end

for i_free = 1:noItems
    patch = set{idx}.item{i_free}.Attributes.patch;
    midx = set{idx}.item{i_free}.Text;
    nurbs = patches{patch}.nurbs;
    outwardPointingNormals = strcmp(set{idx}.Attributes.normal, 'outward');
    inwardPointingNormals = strcmp(set{idx}.Attributes.normal, 'inward');
    bdryNodes = sum(varCol.noCtrlPtsPatch(1:patch-1)) + extractBdryNodes(nurbs,midx,outwardPointingNormals,inwardPointingNormals);
    nodes(counter:counter+numel(bdryNodes)-1) = bdryNodes;
    at = zeros(2,3,'logical');
    at(midx) = true;
    varCol_dummy.dimension = 1;
    
    varCol_dummy.nurbs = subNURBS({nurbs},'at',at.','outwardPointingNormals',outwardPointingNormals,'inwardPointingNormals',inwardPointingNormals,'useFlipping',false);
    
    sub_nurbs(i_free) = varCol_dummy.nurbs;
    varCol_dummy = generateIGAmesh(convertNURBS(varCol_dummy));
    noElemsXiEta = varCol_dummy.patches{1}.noElems;
    index(eBdry:eBdry+noElemsXiEta-1,:) = varCol_dummy.patches{1}.index + repmat(jEl,noElemsXiEta,1);    
    pIndex(eBdry:eBdry+noElemsXiEta-1) = i_free;
    element(eBdry:eBdry+noElemsXiEta-1,:) = counter-1 + varCol_dummy.patches{1}.element;
    elRange_dummy = varCol_dummy.patches{1}.elRange;
    noEl1 = size(elRange_dummy{1},1);
    noEl2 = size(elRange_dummy{2},1);
    noEl3 = size(patches{patch}.elRange{ceil(midx/2)},1);
    elRange{1}(jEl(1)+1:jEl(1)+noEl1,:) = elRange_dummy{1};
    elRange{2}(jEl(2)+1:jEl(2)+noEl2,:) = elRange_dummy{2};
    jEl = jEl + [noEl1,noEl2];
    knotVecs{i_free} = sub_nurbs{i_free}.knots;
    
    temp = repmat(eBdry:eBdry+noElemsXiEta-1,1,noEl3);
    e = sum(noElemsPatch(1:patch-1)) + 1;    
    elemMap(e:e+sum(noElemsPatch(patch))-1) = temp(:);
    neighPatch = NaN;
    for j = 1:numel(connection)
        if connection{j}.Attributes.slave == patch && connection{j}.Attributes.sidx == midx
            neighPatch = connection{j}.Attributes.master;
            sidx = connection{j}.Attributes.midx;
        elseif connection{j}.Attributes.master == patch && connection{j}.Attributes.midx == midx
            neighPatch = connection{j}.Attributes.slave;
            sidx = connection{j}.Attributes.sidx;
        end
    end
    if ~isnan(neighPatch)
        noEl3 = size(patches{neighPatch}.elRange{ceil(sidx/2)},1);
        temp = repmat(eBdry:eBdry+noElemsXiEta-1,1,noEl3);
        e = sum(noElemsPatch(1:neighPatch-1)) + 1;    
        elemMap(e:e+sum(noElemsPatch(neighPatch))-1) = temp(:);
    end
        
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
if 0 % use slow method
    gluedNodes = varCol.gluedNodes;
    for i = 1:length(gluedNodes)
        indices = (nodes(element(:)) == gluedNodes{i}(1));
        if any(indices)
            parentIdx = element(indices);
            for j = 1:length(gluedNodes{i})
                indices = (nodes(element(:)) == gluedNodes{i}(j));
                element(indices) = parentIdx(1);
            end
        end
    end
else % use fast method
    Eps = 1e7*eps;
    controlPts = varCol.controlPts(nodes,:);
    nodesMapBdry = 1:size(controlPts,1);
    [~, gluedNodes] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(norm2(controlPts)), 'OutputAllIndices', true);
    repeatedNode = zeros(numel(gluedNodes),1);
    for i = 1:numel(gluedNodes)
        repeatedNode(i) = numel(gluedNodes{i}) - 1;
    end
    gluedNodes(repeatedNode == 0) = [];

    counter = 1;
    for i = 1:length(gluedNodes)
        childrenNodes_i = gluedNodes{i}(2:end);
        noChildrenNodes_i = numel(childrenNodes_i);
        nodesMapBdry(childrenNodes_i) = gluedNodes{i}(1);
        counter = counter + noChildrenNodes_i;
    end
    element = nodesMapBdry(element);
end
varColBdry.nodes = nodes;
varColBdry.noElems = noElemsBdry;
varColBdry.element = element;
varColBdry.element2 = element2;
varColBdry.index = index;
varColBdry.knotVecs = knotVecs;
varColBdry.n_en = n_en;
varColBdry.noSurfDofs = noSurfDofs;
varColBdry.elemMap = elemMap;
varColBdry.nurbs = sub_nurbs;
varColBdry.degree = degree;
varColBdry.elRange = elRange;
varColBdry.pIndex = pIndex;

function bdryNodes = extractBdryNodes(nurbs,midx,outwardPointingNormals,inwardPointingNormals)
n_xi = nurbs.number(1);
n_eta = nurbs.number(2);
n_zeta = nurbs.number(3);
noControlPts = n_xi*n_eta*n_zeta;
nodes = reshape(1:noControlPts,n_xi,n_eta,n_zeta);
switch midx
    case 1
        bdryNodes = nodes(1,:,:);
        if outwardPointingNormals
            bdryNodes = permute(bdryNodes,[1,3,2]);
        end
    case 2
        bdryNodes = nodes(end,:,:);
        if inwardPointingNormals
            bdryNodes = permute(bdryNodes,[1,3,2]);
        end
    case 3
        bdryNodes = nodes(:,1,:);
        if inwardPointingNormals
            bdryNodes = permute(bdryNodes,[3,2,1]);
        end
    case 4
        bdryNodes = nodes(:,end,:);
        if outwardPointingNormals
            bdryNodes = permute(bdryNodes,[3,2,1]);
        end
    case 5
        bdryNodes = nodes(:,:,1);
        if outwardPointingNormals
            bdryNodes = permute(bdryNodes,[2,1,3]);
        end
    case 6
        bdryNodes = nodes(:,:,end);
        if inwardPointingNormals
            bdryNodes = permute(bdryNodes,[2,1,3]);
        end
end
bdryNodes = bdryNodes(:);






