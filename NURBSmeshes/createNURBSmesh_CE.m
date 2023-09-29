function task = createNURBSmesh_CE(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;
R = varCol{1}.R;
parm = task.msh.parm;
if varCol{1}.boundaryMethod
    fluid = getCEData('R', R, 'parm', parm);
    fluid = makeUniformNURBSDegree(fluid,degree);
    if parm == 2
        Imap{1} = [R*pi/4,R*0.707028554372548];
    else
        Imap = {};
    end
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    t_fluid = task.misc.r_a-R;
    fluid_eye = getOctSphereData('R', R, 'parm', parm, 't', R, 'scaleWeights', false);
    fluid_layer = getEllipsoidData('C', task.misc.r_a, 'parm', parm, 't', t_fluid);
    if parm == 2
        fluid_layer = repeatNURBSknots(insertKnotsInNURBS(fluid_layer,{0.5,0.5,[]}));
    end
    fluid_layer = explodeNURBS(fluid_layer);
    fluid = uniteNURBS({fluid_eye,fluid_layer});
    
    fluid = makeUniformNURBSDegree(fluid,degree);
    if parm == 2
        Imap{1} = [R*pi/4,R*0.707028554372548];
    else
        Imap{1} = [R*pi/2,task.misc.r_a*pi/2];
    end
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
end
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
varCol{1}.nurbs = fluid;
varCol{1}.L_gamma = 2*R;

task.varCol = varCol;
