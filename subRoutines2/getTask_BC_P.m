

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BC'; % BeTSSi submarine
coreMethod = 'IGA';
% coreMethod = {'IGA'};

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

applyLoad = 'radialPulsation';
f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
% f = 5e2; %[1e2 5e2 1e3];             % Frequency
alpha = (0:0.05:360)*pi/180;
method = 'IE';
formulation = 'BGU';
IEbasis = 'Lagrange';
M = 1:2;
degreeElev = 0:1;
calculateSurfaceError = 1;
calculateVolumeError  = 1;
LpOrder = 2; % For error calculation in calcSurfError()
calculateFarFieldPattern = 1;
computeCondNumber = 1;
plotResultsInParaview = 0;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
BC = 'NBC';
N = 3;
loopParameters = {'M','degreeElev','f'};

collectIntoTasks

coreMethod = 'hp_FEM';
degreeElev = 0;
loopParameters = {'M','f'};

collectIntoTasks

coreMethod = 'linear_FEM';
M = 1:3;
collectIntoTasks
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
M = 1:3;
degreeElev = 0;
loopParameters = {'M','f'};
% collectIntoTasks

% M = 1:3;
% degreeElevArr = 0:1;
% coreMethod = {'hp_FEM','h_FEM','C0_IGA'};
% collectIntoTasks
% 
% coreMethod = {'linear_FEM'};
% degreeElevArr = 0;
% M = 1:4;
% collectIntoTasks
% 
% 
% method = {'BEM'};
% formulation = {'CCBIE', 'GBM'}; %, 'CBM', 'GCBIE', 'GBM'};
% M = 1:2;
% degreeElevArr = 0:1;
% collectIntoTasks
% 
% M = 3;
% formulation = {'GBM'}; %, 'CBM', 'GCBIE', 'GBM'};
% plotResultsInParaview = 1;
% collectIntoTasks
% 
% M = 1:2;
% formulation = {'CCBIE', 'GBM'}; %, 'CBM', 'GCBIE', 'GBM'};
% plotResultsInParaview = 0;
% degreeElevArr = 0:1;
% coreMethod = {'hp_FEM','h_FEM','C0_IGA'};
% collectIntoTasks
% 
% coreMethod = {'linear_FEM'};
% degreeElevArr = 0;
% M = 1:3;
% collectIntoTasks
% 