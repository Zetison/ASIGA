function task = createNURBSmesh_Prism(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

x_0 = [0, 0, 0];
task.iem.x_0 = x_0; % The center of the model
alignWithAxis = 'Zaxis';
switch task.misc.method
    case {'IENSG'}
        task.iem.A_2 = eye(3);
        varCol{1}.alignWithAxis = alignWithAxis;
end
L = varCol{1}.L;
if varCol{1}.boundaryMethod
    fluid = getPrismData('x_0',-L/2,'L',L);
    fluid = makeUniformNURBSDegree(fluid,degree);
    refLength = max(L);
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/refLength,{},0,1:3,1,1);
    
    fluid = subNURBS(fluid);
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    t_fluid = varCol{1}.t_fluid;
    fluid = getPrismShellData('x_0',-(L+2*t_fluid)/2,'L',L,'t',t_fluid);
    task.varCol{1}.nurbs = fluid;
    task = defineDomains(task);
    
    refLength = max(L);
    if strcmp(task.misc.method,'PML')
        t_pml = task.pml.t;
        fluid_pml = getPrismShellData('x_0',-(L+2*t_fluid+2*t_pml)/2,'L',L+2*t_fluid,'t',t_pml,'addPMLinfo', true, 't_extraKnots', t_fluid);
        fluid = [fluid,explodeNURBS(fluid_pml,1:3,true)];
        task.varCol{1}.nurbs = fluid;
        task = addPMLtopology(task);
    else
        task.varCol{1}.nurbs = fluid;
    end
    
    fluid = makeUniformNURBSDegree(fluid,degree);
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/refLength,{},0,1:3,1,1);
end
task.varCol{1}.nurbs = fluid;

task.varCol{1}.L_gamma = refLength;