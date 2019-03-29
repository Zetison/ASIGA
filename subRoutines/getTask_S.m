

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'IE'};
formulation = 'BGU';

f = 1e2;             % Frequency
% f = 3.98e4;             % Frequency

M = 1; %3:6

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = pi;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;
calculateVolumeError = 0;
calculateSurfaceError = 1;

calculateFarFieldPattern = 1;
plot3Dgeometry = 0;
degree = 0;
parm = 1;

loopParameters = {'M','parm','method'};
% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation

method = {'BEM'};
formulation = 'CCBIE';
% formulation = 'CBM';
M = 2; %3:6
% formulation = 'CBM';
% coreMethod = {'linear_FEM'};
% formulation = {'CCBIE', 'CBM', 'CHBIE', 'GCBIE', 'GBM', 'GHBIE'};

calculateSurfaceError = 1;
extraGPBEM = 0; % extra quadrature points around singularities for BEM formulations
extraGP = 0; % extra quadrature points around singularities for BEM formulations
% calculateVolumeError = 0;

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = 'BA';
coreMethod = 'IGA';
formulation = 'SL2E';
M = 2;
plot3Dgeometry = 0;
degree = 0:10;
calculateSurfaceError = 1;
parm = 2;
loopParameters = {'degreeElev'};
% computeCondNumber = false;
% loopParameters = {'M','coreMethod','degreeElev','method','formulation'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = '';
M = 5;
plot3Dgeometry = 0;
degree = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
% parm = linspace(0.2,0.6,40);
loopParameters = {'parm','M','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
method = {'RT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degree = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 3:6;
r = 2;

loopParameters = {'parm','M','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
formulation = '';
parm = [];
M = 2:3;
plot3Dgeometry = 0;
degree = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
% collectIntoTasks


