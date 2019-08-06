function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Model2(varCol,parms, M, degree)

solid = NaN;
fluid_i = NaN;


names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

x_0 = [0, 0, 0]; % The origin of the model
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
    fluid = getModel2Data(R_o,t,L,theta2);

    varCol.patchTop = getPatchTopology(fluid);
%     varCol.patchTop{1} = [ones(4,1),zeros(4,1)];
%     varCol.patchTop{1}(2,2) = NaN;
%     varCol.patchTop{1}(4,2) = NaN;
end
L_gamma = L + 2*R_o;

% varCol.chimin = chimin;
% varCol.chimax = chimax;
varCol.L_gamma = L_gamma;
% varCol.Upsilon = Upsilon;