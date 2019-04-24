function [varCol, fluid, solid, fluid_i] = createNURBSmesh_PH(varCol,parms, M, degreeElev)

solid = NaN;
fluid_i = NaN;


names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [-L/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol.method
    case {'IE','IENSG'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        A_2 = [0 1 0;
               0 0 1;
               1 0 0];
        varCol.x_0 = x_0;
        varCol.A_2 = A_2;
        varCol.alignWithAxis = alignWithAxis;
end

if varCol.boundaryMethod
    x_1 = R_2-sqrt(R_2^2-R_1^2);
    c_z = (L+2*x_1)/2;
    c_xy = (R_1+R_2)/2; % 2.5, 3.75; 
    varCol.c_z = c_z;
    varCol.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 22;
    chimax = 23;

    eta1 = R_1*pi/2/(L+(R_1+R_1)*pi/2);
    eta2 = 1-eta1;

    solid = getPHData(R_1, R_2, t, L);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);

    nn = max(2^(M-1)-1,0);

    if M > 0
        solid = insertKnotsInNURBS(solid,{[] linspace2(eta1, eta2, 4) []});
    end

    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});

    fluid = extractOuterSurface(solid);
else
    R_max = max(R_1,R_2);
    s = 0.25 + 0.05*(L-R_max)/R_max;
    c_z = (L+R_1+R_2)/2 + s*R_max;
    c_xy = R_max + s*R_max;

    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

    chimin = 22;
    chimax = 23;

    eta1 = 0.2;
    eta2 = 1-eta1;

    solid = getPHData(R_1, R_2, t, L);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = insertKnotsInNURBS(solid,{[] linspace2(solid.knots{2}(4), solid.knots{2}(6), 4) []});

    nn = max(2^(M-1)-1,0);

    if M > 0
        fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN,[eta1,eta1,linspace2(eta1,eta2,4),eta2,eta2]);
%         fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN);
    end
    fluid = elevateNURBSdegree(fluid,[0 0 1]);

%         fluid = insertKnotsInNURBS(fluid,{[] linspace2(eta1, eta2, 4) []});
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, nn) ...
                                      insertUniform2(fluid.knots{2}, nn) ...
                                      insertUniform2(fluid.knots{3}, nn)});

    fluid = elevateNURBSdegree(fluid,[1 1 1]*degreeElev);

    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});
    solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);
end
L_gamma = L + R_1 + R_2;

varCol.chimin = chimin;
varCol.chimax = chimax;
varCol.L_gamma = L_gamma;
varCol.Upsilon = Upsilon;