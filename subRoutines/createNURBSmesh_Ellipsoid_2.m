
alignWithAxis = 'Zaxis';

% c_z = 26; % 30
% c_xy = 9; % 12,
% c_z = 2;
% c_z = 0.99*(L+2*R_o)/2; % 30
% c_xy = 0.99*R_o; % 2.5, 3.75
switch model
    case 'EL'
        c_z = 6*R_o; % 6*R_o
        c_xy = R_o; % 2.5, 3.75
    case {'SS', 'PS'}
        c_z = R_o; % 30
        c_xy = R_o; % 2.5, 3.75
end
% c_z = c_xy;

x_0 = [0; 0; 0]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)

if 1
    Upsilon = sqrt(c_z^2-c_xy^2);
    chimin = 0.99999999*c_z;
    chimax = 1.00000001*c_z;
else
    Upsilon = 0;
    chimin = 0.99999999*c_xy;
    chimax = 1.00000001*c_z;
end

t = c_z/10;

f_arc = @(s) sqrt(c_xy^2*sin(s).^2+c_z^2*cos(s).^2);
newXiKnots = 2^(i_M-1)-1;
newEtaKnots = round(integral(f_arc,0,pi/2)/(c_xy*pi/2)*(2^(i_M-1)-1));
newZetaKnots = 2^(i_M-1);
% newXiKnots = 0;
% newEtaKnots = 0;
% newZetaKnots = 0;


solid = getEllipsoidalShellData(c_xy,c_xy,c_z,t,alignWithAxis);

solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);
solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, newXiKnots) ...
                                  insertUniform2(solid.knots{2}, newEtaKnots) ...
                                  insertUniform2(solid.knots{3}, newZetaKnots)});


                              
principalLengthXiDir = 2*pi*c_xy;
principalLengthEtaDir = 2*pi*c_z;

L_gamma = 2*c_z;
