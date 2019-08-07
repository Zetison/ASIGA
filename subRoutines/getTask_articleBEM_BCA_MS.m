

scatteringCase = {'MS'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';
 
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = 1e2;             % Frequency
alpha = (0:0.05:180)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
% BC = 'SSBC';
BC = 'SHBC';

method = {'BEM'};
% coreMethod = 'XI';
formulation = {'CCBIE'};
% formulation = {'CRCBIE'};
M = 1:3;
degree = [2,6];
storeSolution = 0;
storeFullVarCol = 0;
agpBEM = 0.6;
solveForPtot = true;
loopParameters = {'method','formulation','M','degree','f','formulation','scatteringCase'};
% collectIntoTasks


f = 1e3;             % Frequency
M = 3;
degree = [5,6];
formulation = {'CBM'};
collectIntoTasks



