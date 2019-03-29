function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model4(varCol, parms, M, degreeElev)

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
    eta1 = R_o/(2*R_o+t*(sqrt(2)-1));
    eta2 = (R_o+sqrt(2)*t)/(2*R_o+t*(sqrt(2)-1));

    principalLengthXiDir = R_o;
    principalLengthEtaDir = R_o;

    L_gamma = R_o;
    fluid = getModel4Data(R_o, t, eta1, eta2);
    fluid = elevateNURBSdegree(fluid,[0 1]);
    fluid = elevateNURBSdegree(fluid,[1 1]*degreeElev);
    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, noNewXiKnots) ...
                                      [insertUniform2([0 eta1], noNewEtaKnots);
                                       insertUniform2([eta1 eta2], 1);
                                       insertUniform2([eta2 1], noNewEtaKnots)]});
end