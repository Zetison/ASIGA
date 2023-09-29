function task = createNURBSmesh_FreeCADsphere(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        error('Not implemented')
end
if varCol{1}.boundaryMethod
    fluid = read_g2('~/kode/ASIGA/miscellaneous/FreeCADsphere.g2');
    surfObjs = [];
    for i = 1:numel(fluid)
        if fluid{i}.d_p == 2
            surfObjs = [surfObjs, i];
        end
    end
    fluid = fluid(surfObjs);
%     fluid = scaleNURBS(fluid,1/1000); % Scale from mm to m
    refLength = 1;
    fluid = makeUniformNURBSDegree(fluid,degree);
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/refLength,{},0);
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    error('Not implemented')
end
varCol{1}.nurbs = fluid;

X = computeBoundingBox(fluid);
varCol{1}.L_gamma = max(X(:,2)-X(:,1));

task.varCol = varCol;