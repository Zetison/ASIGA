function task = createNURBSmesh_Shirron(task)
varCol = task.varCol;
M = task.msh.M;
degree = task.msh.degree;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [-L/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case {'IE','IENSG'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        A_2 = [0 1 0;
               0 0 1;
               1 0 0];
        varCol{1}.x_0 = x_0;
        varCol{1}.A_2 = A_2;
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    c_z = (L+2*R_o)/2;
    c_xy = R_o;
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = c_z;
    chimax = 6.22;

    eta1 = R_o*pi/2/(L+R_o*pi);
    eta2 = 1-eta1;

    solid = getModel3Data(R_o, R_o, t, L, 4);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));

    nn = 2^(M-1)-1;

    solid = insertKnotsInNURBS(solid,{[] linspace2(eta1, eta2, mult-1) []});

    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});

    fluid = extractOuterSurface(solid);

    varCol{1}.patchTop{1} = [ones(4,1),zeros(4,1)];
    varCol{1}.patchTop{1}(2,2) = NaN;
    varCol{1}.patchTop{1}(4,2) = NaN;
else
%     s = 0.25 + 0.05*(L-R_o)/R_o;
%     c_z = (L+2*R_o)/2 + s*R_o;
%     c_xy = R_o + s*R_o;
    c_z = 1.1722*(L+2*R_o)/2; % 30
    c_xy = 2.1212*R_o; % 2.5, 3.75

    Upsilon = sqrt(c_z^2-c_xy^2);
    varCol{1}.r_a = evaluateProlateCoords(0,0,c_z,Upsilon);

%     eta1Arr = [0.35, 0.3, 0.26, 0.25, 0.23, 0.22, 0.22, 0.21, 0.20, 0.19];
%     eta1 = eta1Arr(mult);
%     eta2 = 1-eta1;
    eta1 = 0.21;
    eta2 = 1-eta1;

    chimin = c_z;
    chimax = 6.2;

    solid = getModel3Data(R_o, R_o, t, L, 4);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = insertKnotsInNURBS(solid,{[] linspace2(solid.knots{2}(4), solid.knots{2}(6), mult-1) []});

%         noNewZetaKnots = max(2^(M-1)/2-1,0);
    noNewZetaKnots = max(2^(M-1)/2-1,0);
    if mult < 5
        nn = 2^(M-1)-1;
    else
        nn = 2^(M-1)-1;
    end

    fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN,[eta1,eta1,linspace2(eta1,eta2,mult-1),eta2,eta2]);
%         fluid = embedSolidInEllipsoid(solid,c_z,c_xy,c_xy,alignWithAxis,x_0, 0, NaN);
    fluid = elevateNURBSdegree(fluid,[0 0 1]);

%         fluid = insertKnotsInNURBS(fluid,{[] linspace2(eta1, eta2, mult-1) []});
    fluid = insertKnotsInNURBS(fluid,{insertUniform2(fluid.knots{1}, nn) ...
                                      insertUniform2(fluid.knots{2}, nn) ...
                                      insertUniform2(fluid.knots{3}, noNewZetaKnots)});

    fluid = elevateNURBSdegree(fluid,[1 1 1]*(degree-2));

    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));
end
L_gamma = L + 2*R_o;


varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L_gamma;
varCol{1}.Upsilon = Upsilon;
task.varCol = varCol;