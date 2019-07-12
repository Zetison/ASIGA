function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model4(varCol, parms, M, degree)

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
    eta1 = (R_o-t/2)/(2*R_o+t*sqrt(2));
    eta2 = (R_o-t/2+t*sqrt(2))/(2*R_o+t*sqrt(2));

    principalLengthXiDir = R_o;
    principalLengthEtaDir = R_o;

    L_gamma = R_o;
    if varCol.parm(1) == 1
        noNewXiKnots = 2^M-1; % 8*i_mesh
        noNewEtaKnots = noNewXiKnots;
        noNewEtaKnots2 = max(round(2^(M+1)*t/(R_o*pi/2)-1),0);
        fluid = getModel4Data(R_o, t, eta1, eta2);
        fluid = elevateDegreeInPatches(fluid,[0 1]);
        fluid = elevateDegreeInPatches(fluid,[1 1]*(degree-2));
        fluid{1} = insertKnotsInNURBS(fluid{1},{insertUniform2(fluid{1}.knots{1}, noNewXiKnots) ...
                                          [insertUniform2([0 eta1], noNewEtaKnots);
                                           insertUniform2([eta1 eta2], noNewEtaKnots2);
                                           insertUniform2([eta2 1], noNewEtaKnots)]});
    else
        noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
        noNewEtaKnots = noNewXiKnots;
        noNewEtaKnots2 = max(round(2^(M-1)*t/(R_o*pi/4)-1),0);
        fluid = getModel4Data2(R_o, t);
        fluid = elevateDegreeInPatches(fluid,[1 1]*(degree-2));
        fluid([1:9,16:end]) = insertKnotsInPatches(fluid([1:9,16:end]),noNewXiKnots,noNewEtaKnots);
        fluid(10:15) = insertKnotsInPatches(fluid(10:15),noNewXiKnots,noNewEtaKnots2);
    end
    varCol.patchTop = getPatchTopology(fluid);
end