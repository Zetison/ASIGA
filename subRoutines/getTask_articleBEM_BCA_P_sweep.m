

scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';



plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
% f = [1e2,1e3]; %[1e2 5e2 1e3];             % Frequency
f = linspace(1e2,1e3,1000); %[1e2 5e2 1e3];             % Frequency
% f = 1e2; %[1e2 5e2 1e3];             % Frequency
alpha = beta_s;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
calculateSurfaceError = true;
calculateFarFieldPattern = 0;

BC = 'NBC';

applyLoad = 'radialPulsation';
method = {'BEM'};
formulation = 'CCBIE';
M = 2;
storeSolution = true;
storeFullVarCol = true;
solveForPtot = false;
loopParameters = {'method','M','degree'};
degree = 2;
collectIntoTasks

method = {'BA'};
formulation = 'SL2E';
collectIntoTasks

