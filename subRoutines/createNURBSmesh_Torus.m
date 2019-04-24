function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Torus(varCol, parms, M, degree)

solid = NaN;
fluid_i = NaN;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end
% 
% x_0 = [0, 0, 0];
% switch varCol.method
%     case {'IE','IENSG','ABC'}
%         varCol.x_0 = x_0; % The origin of the model
%         varCol.A_2 = [1 0 0;
%                       0 1 0;
%                       0 0 1];
% end

if varCol.boundaryMethod
    L_gamma = 2*(r_o+r_i);
    solid = getTorusData(r_o,r_i);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*(degree-2));
%     solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, 1) ...
%                                       [] ...
%                                       []});
    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnots = noNewXiKnots;
    noNewZetaKnots = noNewXiKnots;
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                      insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(solid.knots{3}, noNewZetaKnots)});
    fluid = extractOuterSurface(solid);
%     fluid = explodeNURBS(fluid,'eta');
%     fluid = explodeNURBS(fluid,'xi');
    
    varCol.patchTop{1} = [ones(4,1),zeros(4,1)];
end