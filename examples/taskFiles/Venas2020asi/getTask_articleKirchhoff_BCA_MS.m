

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'BCA'; % BeTSSi submarine
misc.coreMethod = 'IGA';
     
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = 1e3;             % Frequency
alpha = (0:0.05:180)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
% BC = 'SSBC';
BC = 'SHBC';

storeSolution = 0;
storeFullVarCol = 0;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
formulation = '';
prePlot.plot3Dgeometry = 0;
degree = 6;
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'misc.method','M','misc.coreMethod'};
% misc.coreMethod = {'linear_FEM'};
misc.coreMethod = {'IGA'};
M = [1,2,3];
collectIntoTasks


