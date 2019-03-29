

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'IE'};
formulation = 'BGU';

f = 20e3;             % Frequency
% f = 3.98e4;             % Frequency

M = 4; %3:6

% alpha_s = 240*pi/180;
% beta_s = 30*pi/180;
% alpha_s = pi;
alpha_s = 0;
beta_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;

loopParameters = {'M','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = '';
M = 5;
plot3Dgeometry = 0;
degreeElev = 0;
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
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 3:5;
r = 2;

loopParameters = {'parm','M','method'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
formulation = '';
parm = [];
M = 2:3;
plot3Dgeometry = 0;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
% collectIntoTasks


