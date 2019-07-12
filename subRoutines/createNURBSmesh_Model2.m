function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model2(varCol,parms, M, degree)

solid = NaN;
fluid_i = NaN;


names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [-L/2-(R_o1-R_o2)/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol.method
    case {'IE','IENSG','MFS'}
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
    c_z = (L+R_o1+R_o2)/2;
%     c_xy = (R_o1+R_o2)/2; % 2.5, 3.75; 
    c_xy = 3*sqrt(665)*(1/20); % 2.5, 3.75; 
    
    varCol.c_z = c_z;
    varCol.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 24.4;
    chimax = 25.7;

    totLength = R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2)+ R_o2*pi/2;
    eta1 = R_o1*pi/2/totLength;
    eta2 = (R_o1*pi/2 + sqrt(L^2+(R_o1-R_o2)^2))/totLength;

    solid = getModel3Data(R_o1, R_o2, t, L);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));

    nn = 2^(M-1)-1;

    solid = insertKnotsInNURBS(solid,{[] linspace2(eta1, eta2, 5) []});

    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, nn) ...
                                      insertUniform2(solid.knots{2}, nn) []});

    fluid = extractOuterSurface(solid);
    varCol.patchTop = getPatchTopology(fluid);
    varCol.patchTop{1} = [ones(4,1),zeros(4,1)];
    varCol.patchTop{1}(2,2) = NaN;
    varCol.patchTop{1}(4,2) = NaN;
end
L_gamma = L + 2*R_o;

varCol.chimin = chimin;
varCol.chimax = chimax;
varCol.L_gamma = L_gamma;
varCol.Upsilon = Upsilon;