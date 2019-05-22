function [varCol, fluid, solid, fluid_i] = createNURBSmesh_Cube(varCol, parms, M, degree)

solid = NaN;
fluid_i = NaN;

names = fieldnames(parms);
for j = 1:numel(names)
    name = names{j};
    eval([name, ' = parms.(names{j});']);
end

varCol.x_0 = [0, 0, 0];
% alignWithAxis = 'Zaxis';
% switch varCol.method
%     case {'IE','IENSG','ABC'}
%         varCol.x_0 = x_0; % The origin of the model
%         varCol.A_2 = [1 0 0;
%                       0 1 0;
%                       0 0 1];
% end

initMeshFactXi = varCol.initMeshFactXi;
if varCol.boundaryMethod
    noNewXiKnots = initMeshFactXi*2^(M-1)-1;
    noNewEtaKnots = initMeshFactXi*2^(M-1)-1;
    
    fluid = getCubeData(a); 
    varCol.patchTop = getPatchTopology(fluid);
%     solid = elevateDegreeInPatches(solid,[0 0 3]);
%     varCol.patchTop = cell(6,1);
%     varCol.patchTop{1} = [2 0;
%                           6 0;
%                           4 0;
%                           5 0];
%     varCol.patchTop{2} = [3 0;
%                           6 3;
%                           1 0;
%                           5 1];
%     varCol.patchTop{3} = [4 0;
%                           6 2;
%                           2 0;
%                           5 2];
%     varCol.patchTop{4} = [1 0;
%                           6 1;
%                           3 0;
%                           5 3];
%     varCol.patchTop{5} = [2 3;
%                           1 0;
%                           4 1;
%                           3 2];
%     varCol.patchTop{6} = [2 1;
%                           3 2;
%                           4 3;
%                           1 0];
    nurbsDegree = fluid{1}.degree(1); % assume all degrees are equal
    
%     solid = explodeNURBS(solid,'eta');
%     solid = explodeNURBS(solid,'xi');
    degree = max(degree,nurbsDegree);
    fluid = elevateDegreeInPatches(fluid,[1 1 1]*(degree-nurbsDegree));
    fluid = insertKnotsInPatches(fluid,noNewXiKnots,noNewEtaKnots);
    L_gamma = a;
end

varCol.L_gamma = L_gamma;