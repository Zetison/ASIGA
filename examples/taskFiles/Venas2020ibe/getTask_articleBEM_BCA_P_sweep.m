

misc.scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'BCA_P'; % BeTSSi submarine
misc.coreMethod = 'IGA';



prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
% f = [1e2,1e3]; %[1e2 5e2 1e3];             % Frequency
f = linspace(1e2,1e3,1000); %[1e2 5e2 1e3];             % Frequency
% f = 1e2; %[1e2 5e2 1e3];             % Frequency
alpha = beta_s;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
err.calculateSurfaceError = true;
calculateFarFieldPattern = 0;

BC = 'NBC';

misc.applyLoad = 'pointPulsation';
misc.method = {'BEM'};
formulation = 'CCBIE';
M = 2;
storeSolution = true;
storeFullVarCol = true;
solveForPtot = false;
loopParameters = {'misc.method','M','degree'};
degree = 2;
collectIntoTasks

misc.method = {'BA'};
formulation = 'SL2E';
collectIntoTasks

