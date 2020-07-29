function varCol = createNURBSmesh_M4(varCol, M, degree)

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
if varCol{1}.boundaryMethod
    L_gamma = R;
    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    noExtraKnots = max(round(2^(M+1)*t/(R*pi/2)-1),0);
    solid = getBeTSSiM4Data('R', R, 't', t, 'parm', varCol{1}.parm);
    solid = makeUniformNURBSDegree(solid,degree);
    fluid = getFreeNURBS(solid);
    
    fluid = insertKnotsInNURBS(fluid,[noNewXiKnots,noNewEtaKnots]);
%     fluid = insertKnotsInNURBS(fluid,[noExtraKnots,noExtraKnots]);
    varCol{1}.patchTop = getPatchTopology(fluid);
end
varCol{1}.nurbs = fluid;
if varCol{1}.useSolidDomain
    varCol{2}.nurbs = solid;
end
if varCol{1}.useInnerFluidDomain
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.L_gamma = L_gamma;