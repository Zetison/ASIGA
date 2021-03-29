function task = createNURBSmesh_M4(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

R = varCol{1}.R;
t = varCol{1}.t;
parm = task.msh.parm;
if varCol{1}.boundaryMethod
    L_gamma = R;
    if parm == 1
    	noExtraKnots = (2^(M-1)-1)/(R*pi/2);
    else
    	noExtraKnots = (2^(M-1)-1)/(R*pi/4);
    end
    solid = getBeTSSiM4Data('R', R, 't', t, 'parm', parm);
    solid = makeUniformNURBSDegree(solid,degree);
    Imap{1} = [R*0.785398163397448, R*0.707028554372548];
    solid = refineNURBSevenly(solid,noExtraKnots,Imap);
    fluid = getFreeNURBS(solid);
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end

varCol{1}.L_gamma = L_gamma;
task.varCol = varCol;