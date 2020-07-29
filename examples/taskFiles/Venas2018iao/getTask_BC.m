

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BC'; % BeTSSi submarine
coreMethod = 'IGA';
% coreMethod = {'IGA'};

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
% f = 5e2; %[1e2 5e2 1e3];             % Frequency
alpha = (0:0.05:360)*pi/180;
method = 'IE';
formulation = 'BGU';
IEbasis = 'Chebyshev';
M = 1:2;
degreeElev = 0;
plotResultsInParaview = 1;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
% BC = 'SSBC';
BC = 'SHBC';
N = 6;
loopParameters = {'M','degreeElev','f'};
% loopParameters = {'f'};

% collectIntoTasks

coreMethod = 'hp_FEM';
degreeElev = 0;
loopParameters = {'M','f'};

% collectIntoTasks

coreMethod = 'linear_FEM';
M = 1:3;
% collectIntoTasks
% N = 4;
% degreeElev = 0;
% M = 1;
% collectIntoTasks

% method = 'IENSG';
% formulation = 'PGU';
% M = 1:4;
% degreeElev = 0:1;
% plotResultsInParaview = 0;
% BC = 'SHBC';
% collectIntoTasks

method = 'BEM';
coreMethod = 'IGA';
formulation = 'CCBIE';
M = 2:3;
degreeElev = 0;
loopParameters = {'M','f'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
% method = {'MFS'};
% f = 5e2;
% M = 1:2;
% degreeElev = 0;
% calculateSurfaceError = 0;
% plotResultsInParaview = 0;
% computeCondNumber = false;
% loopParameters = {'M','method'};
% collectIntoTasks

