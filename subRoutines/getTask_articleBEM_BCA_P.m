

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = [1e2,1e3]; %[1e2 5e2 1e3];             % Frequency
f = 1e3; %[1e2 5e2 1e3];             % Frequency
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
calculateSurfaceError = true;
calculateFarFieldPattern = true;

BC = 'SHBC';

applyLoad = 'radialPulsation';
method = {'BEM'};
formulation = 'CCBIE';
M = 1:2;
storeSolution = true;
storeFullVarCol = true;
loopParameters = {'method','M','degree','f'};
degree = 2:5;
collectIntoTasks

method = {'BA'};
formulation = 'SL2E';
% collectIntoTasks


M = 1;
degree = 3;
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 1;	% Create 30 paraview files in order to visualize a dynamic result
calculateFarFieldPattern = 0;
collectIntoTasks

method = {'BEM'};
formulation = 'CCBIE';
% collectIntoTasks




f = 3e2; %[1e2 5e2 1e3];             % Frequency
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = true;
calculateFarFieldPattern = 0;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = true;
applyLoad = 'radialPulsation';
model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';
method = {'BA'};
formulation = {'SL2E'};
M = 1;
storeSolution = true;
storeFullVarCol = true;
degree = 2;
% degreeElev = 0:8;
loopParameters = {'method','degreeElev'};
% collectIntoTasks

