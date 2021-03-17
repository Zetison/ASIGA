function [nodes1,nodes2,element1,index1,noElems1,pIndex,nodes1_2,nodes2_2] = createSurfaceMesh(varCol1,varCol2)
error('Depricated. Use meshBoundary instead')
patches = varCol1.patches;
noPatches = varCol1.noPatches;
noCtrlPtsPatchSolid = varCol1.noCtrlPtsPatch;
noCtrlPtsPatchFluid = varCol2.noCtrlPtsPatch;
gluedNodes1 = varCol1.gluedNodes;
gluedNodes2 = varCol2.gluedNodes;
e = 1;
nodes1 = [];
nodes2 = [];
index1 = [];
noElems1 = 0;
element1 = [];
noCtrlPtsSurf = 0;
outerBoundary = true;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    nodes1_loc = zeros(1,n_xi*n_eta);
    nodes2_loc = zeros(1,n_xi*n_eta);
    counter = 1;
    for j = 1:n_eta
        for i = 1:n_xi
            nodes1_loc(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
            nodes2_loc(counter) = n_xi*(j-1) + i;
            counter = counter + 1;
        end
    end
    nodes1 = [nodes1, nodes1_loc+sum(noCtrlPtsPatchSolid(1:patch-1))];
    nodes2 = [nodes2, nodes2_loc+sum(noCtrlPtsPatchFluid(1:patch-1))];

    nurbs = patches{patch}.nurbs;
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = subNURBS({nurbs},'at',[0,0;0,0;1-outerBoundary,outerBoundary]);
    varCol_dummy = generateIGAmesh(convertNURBS(varCol_dummy));
    noElems1_loc = varCol_dummy.patches{1}.noElems;
    index1_loc = varCol_dummy.patches{1}.index;   
    element1_loc = varCol_dummy.patches{1}.element;
    
    
    noElems1 = noElems1 + noElems1_loc;
    element1 = [element1; element1_loc+noCtrlPtsSurf];
    index1 = [index1; index1_loc];
    pIndex(e:e+noElems1_loc-1) = patch;
    e = e + noElems1_loc;
    noCtrlPtsSurf = noCtrlPtsSurf + n_xi*n_eta;
end
nodes1_2 = nodes1;
nodes2_2 = nodes2;

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes1)
    parentIdxSolid = gluedNodes1{i}(1);
    for j = 2:length(gluedNodes1{i})
        indices = (nodes1 == gluedNodes1{i}(j));
        nodes1(indices) = parentIdxSolid;
    end
end
for i = 1:length(gluedNodes2)
    parentIdxFluid = gluedNodes2{i}(1);
    for j = 2:length(gluedNodes2{i})
        indices = (nodes2 == gluedNodes2{i}(j));
        nodes2(indices) = parentIdxFluid;
    end
end