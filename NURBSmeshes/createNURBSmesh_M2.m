function varCol = createNURBSmesh_M2(varCol, M, degree)

R = varCol{1}.R;
t = varCol{1}.t;
L = varCol{1}.L;
theta1 = varCol{1}.theta1;
theta2 = varCol{1}.theta2;
if varCol{1}.boundaryMethod
    fluid = getBeTSSiM2Data('R', R, 'L', L, 't', t, 'parm', varCol{1}.parm, 'theta1', theta1, 'theta2', theta2);
    fluid = makeUniformNURBSDegree(fluid,degree);
    Imap{1} = [0.010000018518611 0.000023554074421   0.000033333395062   0.007066209240776   0.009976464444190 0.009999962962798];
    Imap{2} = [4.681336475976538 2.985867581518450   2.989983333287038   t/2*pi/2 0.015641296477826 0.011052467370932];
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2),Imap);
    varCol{1}.patchTop = getPatchTopology(fluid);
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
end
varCol{1}.nurbs = fluid;
varCol{1}.L_gamma = L + 2*R;

