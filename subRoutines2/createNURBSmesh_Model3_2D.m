% L  = 20;
alignWithAxis = 'Zaxis';
degreeElev = 0;
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

xi1 = 0.1;
xi2 = 0.37;
xi4 = 1-xi2;
xi5 = 1-xi1;

i_mesh1 = ceil(pi/2*R_o1*(2^i_mesh-1));
i_mesh2 = ceil(0.5*L*(2^i_mesh-1));
i_mesh3 = ceil(pi/2*R_o2*(2^i_mesh-1));

noNewEtaKnots = 2*(2^i_mesh-1)*3;

% solid = getModel2Data_2D(R_o, L, t, x_0);
solid = getModel3Data_2D(R_o2, R_o1, L, x_0);
nurbs = embedSolidInEllipse(solid,c_x,c_y);

nurbs = elevateNURBSdegree(nurbs,[0 1]);

nurbs = insertKnotsInNURBS(nurbs,{[linspace2(0,xi1,i_mesh1) linspace2(xi1,0.25,i_mesh2) linspace2(0.25,xi2,i_mesh2) linspace2(xi2,0.5,i_mesh3) ...
                                   linspace2(0.5,xi4,i_mesh3) linspace2(xi4,0.75,i_mesh2) linspace2(0.75,xi5,i_mesh2) linspace2(xi5,1,i_mesh1)] ...
                                   insertUniform2(nurbs.knots{2}, noNewEtaKnots)});

fluid_o = nurbs;

principalLengthXiDir = 2*L+2*R_o;