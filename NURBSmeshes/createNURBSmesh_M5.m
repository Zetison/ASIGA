function varCol = createNURBSmesh_M5(varCol, M, degree)

x_0 = [0, 0, 0];
switch varCol{1}.method
    case {'IE','IENSG','ABC'}
        varCol{1}.x_0 = x_0; % The origin of the model
        varCol{1}.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
        varCol{1}.alignWithAxis = alignWithAxis;
end


if varCol{1}.boundaryMethod
    principalLengthXiDir = R_o;
    principalLengthEtaDir = R_o;

    L_gamma = R_o;
    eta2 = (R_o+L)/(L+2*R_o);
    eta1 = R_o/(L+2*R_o);
    switch varCol{1}.model
        case 'M5A'
            fluid = getModel5Data(R_o, eta1, eta2, L, l, 'Xaxis');
        case 'M5B'
            fluid = getModel5Data(R_o, eta1, eta2, L, l, 'Zaxis');
    end

    fluid = elevateDegreeInPatches(fluid,[1 1]*(degree-2));
    newKnots = {[] linspace2(eta1,eta2,round(L/R_o))};
    fluid{1} = insertKnotsInNURBS(fluid{1},newKnots);
    fluid{2} = insertKnotsInNURBS(fluid{2},newKnots);
    
    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
    varCol{1}.patchTop = getPatchTopology(fluid);
end
varCol{1}.nurbs = fluid;
if varCol{1}.useSolidDomain
    varCol{2}.nurbs = solid;
end
if varCol{1}.useInnerFluidDomain
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
varCol{1}.Upsilon = Upsilon;