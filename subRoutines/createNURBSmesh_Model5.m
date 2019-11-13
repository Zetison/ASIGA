function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model5(varCol, parms, M, degree)

solid = NaN;
fluid_i = NaN;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [0, 0, 0];
switch varCol.method
    case {'IE','IENSG','ABC'}
        varCol.x_0 = x_0; % The origin of the model
        varCol.A_2 = [1 0 0;
                      0 1 0;
                      0 0 1];
end


if varCol.boundaryMethod
    principalLengthXiDir = R_o;
    principalLengthEtaDir = R_o;

    L_gamma = R_o;
    eta2 = (R_o+L)/(L+2*R_o);
    eta1 = R_o/(L+2*R_o);
    switch model
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
    varCol.patchTop = getPatchTopology(fluid);
end