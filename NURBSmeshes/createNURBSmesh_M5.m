function varCol = createNURBSmesh_M5(varCol, M, degree)

R = varCol{1}.R;
l = varCol{1}.l;
L = varCol{1}.L;
if varCol{1}.boundaryMethod    
    fluid = getBeTSSiM5Data('R', R, 'type', varCol{1}.type, 'L', L, 'l', l, 'parm', varCol{1}.parm);

    fluid = makeUniformNURBSDegree(fluid,degree);
    fluid = refineNURBSevenly(fluid,2.2*10/3);
    fluid = insertKnotsInNURBS(fluid,(2^(M-1)-1)*[1,1]);
    
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
end
varCol{1}.nurbs = fluid;

varCol{1}.L_gamma = L;