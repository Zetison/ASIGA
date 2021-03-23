function task = createNURBSmesh_Cube(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;


task.iem.x_0 = [0, 0, 0];

initMeshFactXi = varCol{1}.initMeshFactXi;
if varCol{1}.boundaryMethod
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = initMeshFactXi*2^(M-1)-1;
    
    fluid = getCubeData(a); 
    varCol{1}.patchTop = getPatchTopology(fluid);
    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    
    degree = max(degree,nurbsDegree);
    fluid = elevateDegreeInPatches(fluid,[1 1 1]*(degree-nurbsDegree));
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
    L_gamma = a;
end

varCol{1}.nurbs = fluid;

varCol{1}.L_gamma = L_gamma;
task.varCol = varCol;