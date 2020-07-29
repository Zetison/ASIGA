scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = [1e2 1e3];
f = 1e2;
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
BC = 'SHBC';

method = {'BEM'};
formulation = {'CCBIE'};
% formulation = {'CBM'};
M = 1:3;
% M = 1;
degree = [2,6];
% degree = 2;
storeSolution = 0;
storeFullVarCol = 0;
% agpBEM = 0.6;
agpBEM = 1.4;
solveForPtot = true;
loopParameters = {'method','formulation','M','degree','f','scatteringCase'};

collectIntoTasks

f = 1e3;
degree = 6;
M = 3;
formulation = {'CBM'};
collectIntoTasks


formulation = {'GCBIE'};
extraGPBEM = 200; % extra quadrature points around singularities for BEM formulations
M = 1;
% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calculateSurfaceError = 1;
applyLoad = 'pointPulsation';
BC = 'NBC';
M = 1:3;
degree = [2,5];
method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks

