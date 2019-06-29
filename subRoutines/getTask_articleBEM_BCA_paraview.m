scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

alpha = (0:0.05:360)*pi/180;

BC = 'SHBC';

method = {'BEM'};
formulation = 'CCBIEC';


f = 1e3;             % Frequency
M = 3;
degree = 5;
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = 0;
calculateFarFieldPattern = 0;
plotTimeOscillation   = 1;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = 0;
storeSolution = 0;
storeFullVarCol = 0;
agpBEM = 0.6;
solveForPtot = true;
loopParameters = {'method','degree'};
collectIntoTasks

