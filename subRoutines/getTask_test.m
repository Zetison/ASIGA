

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'IE'};
formulation = 'BGU';
% k = 1;
% c = 1524;
% omega = c*k;
% f = omega/(2*pi);
f = 1e3;             % Frequency
% f = 3.98e4;             % Frequency

M = 2; %3:6
N = 6;
IEbasis = 'Chebyshev';

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha_s = 0;
beta_s = pi/2;
% alpha_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;
calculateVolumeError = 0;
calculateSurfaceError = 1;

calculateFarFieldPattern = 1;
plot3Dgeometry = 0;
degree = 4;
parm = 2;

loopParameters = {'M','parm','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = 'SS'; % PS = point solution, SS = spherical solution
% formulation = 'PS'; % PS = point solution, SS = spherical solution
M = 1;

calculateSurfaceError = 1;
computeCondNumber = false;
delta = [0.1,0.2,0.4,0.6,0.8];
extraGP = [1,2,4,8,16,32];
loopParameters = {'delta','extraGP','method'};

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% BA simulation
% method = {'BA'};
% coreMethod = 'IGA';
% formulation = 'SL2E';
% M = 1:5;
% plot3Dgeometry = 0;
% degreeElev = 0;
% calculateSurfaceError = 1;
% parm = [1,2];
% loopParameters = {'M','parm','method'};
% % computeCondNumber = false;
% % loopParameters = {'M','coreMethod','degreeElev','method','formulation'};
% collectIntoTasks
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% MFS simulation
% method = {'MFS'};
% formulation = '';
% M = 5;
% plot3Dgeometry = 0;
% degreeElev = 0;
% calculateSurfaceError = 0;
% computeCondNumber = false;
% % parm = linspace(0.2,0.6,40);
% loopParameters = {'parm','M','method'};
% % collectIntoTasks
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% RT simulation
% method = {'RT'};
% formulation = '';
% M = 3;
% plot3Dgeometry = 0;
% degreeElev = 0;
% calculateSurfaceError = 0;
% computeCondNumber = false;
% plotFarField = 1;
% applyLoad = 'planeWave';
% parm = 3:6;
% r = 2;
% 
% loopParameters = {'parm','M','method'};
% % collectIntoTasks
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% KDT simulation
% method = {'KDT'};
% formulation = '';
% parm = [];
% M = 2:3;
% plot3Dgeometry = 0;
% degreeElev = 0;
% calculateSurfaceError = 0;
% computeCondNumber = false;
% loopParameters = {'M','method'};
% % collectIntoTasks
% 

