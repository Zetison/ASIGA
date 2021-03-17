

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'BCA_P'; % BeTSSi submarine
misc.coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = [1e2 1e3];             % Frequency
% f = 1e3;             % Frequency
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
err.calculateSurfaceError = true;
calculateFarFieldPattern = false;

BC = 'NBC';

misc.applyLoad = 'pointPulsation';
misc.method = {'BA'};
formulation = 'SL2E';
M = 1:7;
% M = 1:5;
storeSolution = 0;
storeFullVarCol = 0;
solveForPtot = true;
loopParameters = {'misc.method','M','degree','f'};
degree = 2:4;
collectIntoTasks


