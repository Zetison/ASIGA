function varCol = createNURBSmesh_M1(varCol, M, degree)

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
if varCol{1}.boundaryMethod
    fluid = getBeTSSiM1Data('R', R, 't', t, 'parm', varCol{1}.parm);
    fluid = makeUniformNURBSDegree(fluid,degree);
    fluid = refineNURBSevenly(fluid,0.5);
    fluid = insertKnotsInNURBS(fluid,(2^(M-1)-1)*[1,1]);
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

