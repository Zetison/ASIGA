function task = createNURBSmesh_M31(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

R1 = varCol{1}.R1;
R2 = varCol{1}.R2;
t = varCol{1}.t;
L = varCol{1}.L;
varCol{1}.x_0 = [-L/2-(R2-R1)/2, 0, 0]; % The origin of the model
if varCol{1}.boundaryMethod

    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', 2);
    solid = makeUniformNURBSDegree(solid,degree);
    Imap{1} = [R2*pi/2,R1*pi/2,(R2-t)*pi/2,(R1-t)*pi/2];
    solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R2*pi/2),Imap);

    fluid = subNURBS(solid, 'at', [0,0;0,0;0,1]);
    varCol{1}.patchTop = getPatchTopology(fluid);
    if numel(varCol) > 2
        fluid_i_outer = subNURBS(solid, 'at', [0,0;0,0;1,0]);

        R = varCol{3}.R;
        t2 = varCol{3}.t;
        L2 = varCol{3}.L;
        
        fluid_i_inner = getBeTSSiM1Data('R', R, 'L', L2, 't', t2);
        fluid_i_inner = translateNURBS(fluid_i_inner,[-(L-L2),0,0]);
        fluid_i_inner = makeUniformNURBSDegree(fluid_i_inner,degree);
        Imap{1} = [R*pi/2 4.242171326235288];
        Imap{2} = [R*pi/4 2.121085663117643];
        fluid_i_inner = refineNURBSevenly(fluid_i_inner,(2^(M-1)-1)/(R*pi/2),Imap);
        fluid_i = [fluid_i_inner, fluid_i_outer];

        varCol{3}.patchTop = getPatchTopology(fluid_i);

        varCol_fluid_i_inner.dimension = 1;
        varCol_fluid_i_inner.nurbs = fluid_i_inner;
        varCol_fluid_i_inner = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i_inner)));

        varCol_fluid_i_outer.dimension = 1;
        varCol_fluid_i_outer.nurbs = fluid_i_outer;
        varCol_fluid_i_outer = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_fluid_i_outer)));

        varCol{3}.elemsOuter = (1:varCol_fluid_i_outer.noElems)+varCol_fluid_i_inner.noElems;
        varCol{3}.noDofsInner = varCol_fluid_i_inner.noDofs;
        varCol{3}.noElemsInner = varCol_fluid_i_inner.noElems;
    end
    
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
end
if degree < 4
    warning(['Using degree=4 instead of degree=' num2str(degree)])
end
varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.L_gamma = L + R1 + R2;
task.varCol = varCol;