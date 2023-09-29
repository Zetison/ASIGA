function task = createNURBSmesh_Torus(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

if varCol.boundaryMethod
    L_gamma = 2*(r_o+r_i);
    solid = getTorusData(r_o,r_i);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));
    
    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    noNewZetaKnots = noNewXiKnots;
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                      insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(solid.knots{3}, noNewZetaKnots)});
    fluid = extractOuterSurface(solid);
%     fluid = explodeNURBS(fluid,'eta');
%     fluid = explodeNURBS(fluid,'xi');
    
    varCol.patchTop = getPatchTopology(fluid);
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end
task.varCol = varCol;