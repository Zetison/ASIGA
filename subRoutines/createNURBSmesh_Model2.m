
alignWithAxis = 'Zaxis';

c_y = 2*R_o; % 30
c_x = L/2+2*R_o; % 12
% c_y = c_x;

x_0 = [L/2; 0]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)


f = sqrt(c_x^2-c_y^2);
[R_a, ~] = evaluateEllipticCoords(c_x,0,f);
rm = zeros(N+2,1);
for m = 1:N+2
    rm(m) = R_a + (m-1)*R_a;
end

xi1 = 0.14;
xi2 = 0.37;
xi3 = 0.4;
xi4 = 0.45;
xi5 = 0.55;
xi6 = 0.6;
xi7 = 0.63;
xi8 = 0.86;

i_mesh1 = ceil(pi/2*R_o*(2^i_mesh-1));
i_mesh2 = ceil(0.5*L*(2^i_mesh-1));
i_mesh3 = ceil(R_o*(2^i_mesh-1));
i_mesh4 = ceil(R_o*(2^i_mesh-1));
i_mesh5 = ceil(0.5*6*t*(2^i_mesh-1));
noNewEtaKnots = 2*(2^i_mesh-1)*3;

solid = getModel2Data_2D(0.5*R_o, R_o, L, t, x_0);
nurbs = embedSolidInEllipse(solid,c_x,c_y);

nurbs = elevateNURBSdegree(nurbs,[0 1]);

nurbs = insertKnotsInNURBS(nurbs,{[linspace2(0,xi1,i_mesh1) linspace2(xi1,0.25,i_mesh2) linspace2(0.25,xi2,i_mesh2) linspace2(xi2,xi3,i_mesh3) ...
                                   linspace2(xi3,xi4,i_mesh4) linspace2(xi4,0.5,i_mesh5) linspace2(0.5,xi5,i_mesh5) linspace2(xi5,xi6,i_mesh4) ...
                                   linspace2(xi6,xi7,i_mesh3) linspace2(xi7,0.75,i_mesh2) linspace2(0.75,xi8,i_mesh2) linspace2(xi8,1,i_mesh1)] ...
                                   insertUniform2(nurbs.knots{2}, noNewEtaKnots)});

fluid_o = nurbs;

principalLengthXiDir = 2*L+2*R_o;