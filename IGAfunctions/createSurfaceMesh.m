function [solidNodes,fluidNodes,solidXiEtaMesh,solidIndexXiEta,solidNoElemsXiEta,pIndex,solidNodes2,fluidNodes2] ...
                        = createSurfaceMesh(varColSolid,varColFluid)
patches = varColSolid.patches;
noPatches = varColSolid.noPatches;
knotVecs = varColSolid.knotVecs;
noCtrlPtsPatchSolid = varColSolid.noCtrlPtsPatch;
noCtrlPtsPatchFluid = varColFluid.noCtrlPtsPatch;
gluedNodesSolid = varColSolid.gluedNodes;
gluedNodesFluid = varColFluid.gluedNodes;
e = 1;
solidNodes = [];
fluidNodes = [];
solidIndexXiEta = [];
solidNoElemsXiEta = 0;
solidXiEtaMesh = [];
noCtrlPtsSurf = 0;
outerBoundary = true;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    p_xi = patches{patch}.nurbs.degree(1);
    p_eta = patches{patch}.nurbs.degree(2);
    solidNodes_loc = zeros(1,n_xi*n_eta);
    fluidNodes_loc = zeros(1,n_xi*n_eta);
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    counter = 1;
    for j = 1:n_eta
        for i = 1:n_xi
            solidNodes_loc(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
            fluidNodes_loc(counter) = n_xi*(j-1) + i;
            counter = counter + 1;
        end
    end
    solidNodes = [solidNodes, solidNodes_loc+sum(noCtrlPtsPatchSolid(1:patch-1))];
    fluidNodes = [fluidNodes, fluidNodes_loc+sum(noCtrlPtsPatchFluid(1:patch-1))];

    nurbs = patches{patch}.nurbs;
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = subNURBS({nurbs},'at',[0,0;0,0;1-outerBoundary,outerBoundary]);
    varCol_dummy = generateIGAmesh(convertNURBS(varCol_dummy));
    solidNoElemsXiEta_loc = varCol_dummy.patches{1}.noElems;
    solidIndexXiEta_loc = varCol_dummy.patches{1}.index;   
    solidXiEtaMesh_loc = varCol_dummy.patches{1}.element;
    
    
    solidNoElemsXiEta = solidNoElemsXiEta + solidNoElemsXiEta_loc;
    solidXiEtaMesh = [solidXiEtaMesh; solidXiEtaMesh_loc+noCtrlPtsSurf];
    solidIndexXiEta = [solidIndexXiEta; solidIndexXiEta_loc];
    pIndex(e:e+solidNoElemsXiEta_loc-1) = patch;
    e = e + solidNoElemsXiEta_loc;
    noCtrlPtsSurf = noCtrlPtsSurf + n_xi*n_eta;
end
solidNodes2 = solidNodes;
fluidNodes2 = fluidNodes;

% Glue nodes in 2D mesh
for i = 1:length(gluedNodesSolid)
    parentIdxSolid = gluedNodesSolid{i}(1);
    for j = 2:length(gluedNodesSolid{i})
        indices = (solidNodes == gluedNodesSolid{i}(j));
        solidNodes(indices) = parentIdxSolid;
    end
end
for i = 1:length(gluedNodesFluid)
    parentIdxFluid = gluedNodesFluid{i}(1);
    for j = 2:length(gluedNodesFluid{i})
        indices = (fluidNodes == gluedNodesFluid{i}(j));
        fluidNodes(indices) = parentIdxFluid;
    end
end