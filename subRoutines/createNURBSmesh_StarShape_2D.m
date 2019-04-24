% L  = 20;
alignWithAxis = 'Zaxis';
degreeElev = 0;
c_y = 1.3; % 30
c_x = 1.3; % 12
% c_y = c_x;

x_0 = [0; 0]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)


f = sqrt(c_x^2-c_y^2);
[R_a, ~] = evaluateEllipticCoords(c_x,0,f);
rm = zeros(N+2,1);
for m = 1:N+2
    rm(m) = R_a + (m-1)*R_a;
end


noNewXiKnots = 2*(2^i_mesh-1)*3;
noNewEtaKnots = 2*(2^i_mesh-1)*3;

% solid = getModel2Data_2D(R_o, L, t, x_0);
solid = getStarShapeData_2D();

solid = elevateNURBSdegree(solid,1);
nurbs = embedSolidInEllipse(solid,c_x,c_y);

nurbs = elevateNURBSdegree(nurbs,[0 1]);

nurbs = insertKnotsInNURBS(nurbs,{insertUniform2(nurbs.knots{1}, noNewXiKnots) ...
                                   insertUniform2(nurbs.knots{2}, noNewEtaKnots)});

fluid_o = nurbs;

principalLengthXiDir = 2*L+2*R_o;