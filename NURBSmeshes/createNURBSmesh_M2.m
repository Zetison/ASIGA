function varCol = createNURBSmesh_M2(varCol, M, degree)

x_0 = [0, 0, 0];
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.x_0 = x_0; % The origin of the model
        varCol{1}.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
        varCol{1}.alignWithAxis = alignWithAxis;
end

R = varCol{1}.R;
t = varCol{1}.t;
L = varCol{1}.L;
theta1 = varCol{1}.theta1;
theta2 = varCol{1}.theta2;
if varCol{1}.boundaryMethod
    fluid = getBeTSSiM2Data('R', R, 't', t, 'parm', varCol{1}.parm, 'theta1', theta1, 'theta2', theta2);
    fluid = makeUniformNURBSDegree(fluid,degree);
    fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R*pi/2));
    varCol{1}.patchTop = getPatchTopology(fluid);
end
varCol{1}.nurbs = fluid;
if varCol{1}.useSolidDomain
    varCol{2}.nurbs = solid;
end
if varCol{1}.useInnerFluidDomain
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.L_gamma = L + 2*R;

