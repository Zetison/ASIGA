

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M4'; % BeTSSi misc.model 5A and BeTSSi misc.model 5B

misc.method = {'BEM'};
formulation = {'CBM','CCBIE'};
% formulation = {'CCBIE'};

f = 10e3;             % Frequency
% f = 3e3;             % Frequency

M = 7;
% M = 1;
parm = [2,1];
parm = 2;
degree = 2;
beta_s = 0;
alpha_s = 30*pi/180;

prePlot.plot3Dgeometry = 0;
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
err.err.calculateSurfaceError = 0;
calculateFarFieldPattern = 1;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = 0;
storeSolution = 0;
storeFullVarCol = 0;
solveForPtot = true;
loopParameters = {'misc.method','M','degree'};
% collectIntoTasks


misc.method = {'KDT'};
solveForPtot = false;
formulation = {''};
collectIntoTasks
