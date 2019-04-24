alignWithAxis = 'Zaxis';

c_z = 0.99*(L+2*R_o)/2; % 30
c_xy = 0.99*R_o; % 2.5, 3.75

eta2 = 0.7;
eta1 = 0.30;

x_0 = [0; 0; L/2]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)


Upsilon = sqrt(c_z^2-c_xy^2);

chimin = 5.93;
chimax = 6.2;


solid = getModel3Data(R_o, R_o, R_i, R_i, L, alignWithAxis, eta1, eta2, x_0);
solid = elevateNURBSdegree(solid,[0 0 1]);
solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);

noNewXiKnots = 2^(i_mesh-1)-1; % 8*i_mesh
noNewEtaKnotsEnds = noNewXiKnots;
noNewEtaKnotsMiddle = floor(noNewXiKnots*L/(R_o*pi/2));

noNewZetaKnots = 2^(i_mesh-1)-1;

solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                  [insertUniform2([0 eta1], noNewEtaKnotsEnds); ...
                                   insertUniform2([eta1 eta2], noNewEtaKnotsMiddle); ...
                                   insertUniform2([eta2 1], noNewEtaKnotsEnds)] ...
                                  insertUniform2(solid.knots{3}, noNewZetaKnots)});



if useInnerFluidDomain
    fluid_i = getModel3InternalWaterData(R_i,R_i,L,alignWithAxis, eta1, eta2, x_0);
    fluid_i = elevateNURBSdegree(fluid_i,[0 0 1]);
    fluid_i = elevateNURBSdegree(fluid_i,[1 1 1]*degreeElev);
    fluid_i = insertKnotsInNURBS(fluid_i,{insertUniform2(fluid_i.knots{1}, noNewXiKnots) ...
                                      [insertUniform2([0 eta1], noNewEtaKnotsEnds); ...
                                       insertUniform2([eta1 eta2], noNewEtaKnotsMiddle); ...
                                       insertUniform2([eta2 1], noNewEtaKnotsEnds)] ...
                                      insertUniform2(fluid_i.knots{3}, noNewZetaKnots)});
end

principalLengthXiDir = 2*pi*R_o;
principalLengthEtaDir = R_o*pi + L;

L_gamma = L + R_o + R_o;