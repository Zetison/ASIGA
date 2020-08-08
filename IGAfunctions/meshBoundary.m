function [nodes, noElems, element, element2, index, pIndex, n_en, noSurfDofs] = meshBoundary(varCol,outerBoundary)
gluedNodes = varCol.gluedNodes;

degree = varCol.degree(1:2);
patches = varCol.patches;
knotVecs = varCol.knotVecs;
noPatches = varCol.noPatches;

p_xi = degree(1); % assume p_xi is equal in all patches
p_eta = degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);


noParams = 2;
noSurfDofs = 0;
noElemsPatch = zeros(noPatches,1);
noEl = zeros(noPatches,noParams);
for i = 1:noPatches
    n_xi = patches{i}.nurbs.number(1);
    n_eta = patches{i}.nurbs.number(2);
    noSurfDofs = noSurfDofs + n_xi*n_eta;
    for j = 1:noParams
        noEl(i,j) = size(patches{i}.elRange{j},1);
    end
    noElemsPatch(i) = size(patches{i}.elRange{1},1)*size(patches{i}.elRange{2},1);
end

nodes = zeros(1,noSurfDofs);
counter = 1;
shiftIdx = 0;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    for j = 1:n_eta
        for i = 1:n_xi
            if outerBoundary
                nodes(counter) = shiftIdx + (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
            else
                nodes(counter) = shiftIdx + n_xi*(j-1) + i;
            end
            counter = counter + 1;
        end
    end
    shiftIdx = shiftIdx + n_zeta*n_eta*n_xi;
end
noElems = sum(noElemsPatch);
pIndex = zeros(noElems,1);
element = zeros(noElems,n_en);
index = zeros(noElems,noParams);
e = 1;
maxDof = 0;
jEl = zeros(1,2);
for i = 1:noPatches
    nurbs = patches{i}.nurbs;
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = subNURBS({nurbs},'at',[0,0;0,0;1-outerBoundary,outerBoundary]);
    varCol_dummy = generateIGAmesh(convertNURBS(varCol_dummy));
    noElemsXiEta = varCol_dummy.patches{1}.noElems;
    index(e:e+noElemsPatch(i)-1,:) = varCol_dummy.patches{1}.index + repmat(jEl,noElemsXiEta,1);    
    pIndex(e:e+noElemsXiEta-1) = i;
    element(e:e+noElemsXiEta-1,:) = maxDof + varCol_dummy.patches{1}.element;
    jEl = jEl + noEl(i,:);
    maxDof = maxDof + n_xi*n_eta;
    e = e + noElemsPatch(i);
end
element2 = element;
% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    indices = (nodes(element(:)) == gluedNodes{i}(1));
    parentIdx = element(indices);
    element(indices) = parentIdx;
    for j = 2:length(gluedNodes{i})
        indices = (nodes(element(:)) == gluedNodes{i}(j));
        element(indices) = parentIdx;
    end
end












