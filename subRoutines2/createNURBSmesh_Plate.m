
x_0 = [0 0 0]; % The origin of the model
Upsilon = 0;

fluid_o = getRectangleIn3DData([0 0 0],[1 0 0], [1 1 0], [0 1 0]);
fluid_o = elevateNURBSdegree(fluid_o,[1 1]*degreeElev);

noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
noNewEtaKnots = noNewXiKnots;

fluid_o = insertKnotsInNURBS(fluid_o,{insertUniform2(fluid_o.knots{1}, noNewXiKnots) ...
                                  insertUniform2(fluid_o.knots{2}, noNewEtaKnots)});

principalLengthXiDir = 1;
principalLengthEtaDir = 1;

L_gamma = 1;