
x_0 = [0; 0; L/2];

x_0wave = x_0;


alignWithAxis = 'Zaxis';
eta2 = 0.765;
eta1 = 0.235;

stretchingFact = 1.30;

c_z = 0.99*L/2; % 30
c_xy = 0.99*R_o; % 12
Upsilon = sqrt(c_z^2 - c_xy^2);
chimin = 9.8;
chimax = 11.2;

solid = getBarrelData(R_o,R_i,eta1,eta2,L,alignWithAxis, x_0);

noNewXiKnots = 2^(i_M-1)-1;  
noNewZetaKnots = 2^(i_M-1)-1;
i_M1 = 2^(i_M-1)-1;
i_M2 = (2^(i_M-1)-1)*L/R_o;
i_M3 = 2^(i_M-1)-1;


solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);

solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                  [insertUniform2([0 eta1], i_M1); ...
                                   insertUniform2([eta1 eta2], i_M2);
                                   insertUniform2([eta2 1], i_M3)] ...
                                  insertUniform2(solid.knots{3}, noNewZetaKnots)});

                              
                       
principalLengthXiDir = 2*pi*R_o;
principalLengthEtaDir = 2*R_o + L;
L_gamma = L;


if ~strcmp(BC, 'SHBC')
    error('This is not implemented')
end
    
    
    