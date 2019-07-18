

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M3'; % BeTSSi model 5A

method = {'IENSG','IE'};
method = {'IENSG'};
method = {'IE'};
formulation = 'BGU';

f = 1e3; %[1e2 5e2 1e3];             % Frequency
omega = 2*pi*f;
k = omega/1500;
parm = 1;
M = 3;
degree = 2;
N = [3,5];
N = 3;
loopParameters = {'degree','M','N','method'};
alpha = (0:0.1:180)*pi/180;
% alpha = (87:0.1:91)*pi/180;
% alpha = sort([alpha, pi/2-atan((5-3)/41), 3*pi/2+atan((5-3)/41)]);
alpha = sort([alpha, pi/2-atan((5-3)/41)]);
% alpha = 88*pi/180;

plot3Dgeometry = 0;

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
% collectIntoTasks

method = {'KDT'};
M = 5;
degree = 2;
loopParameters = {'degree','M','method'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = 'SS'; % PS = point solution, SS = spherical solution
M = 1:3;
parm = linspace(0.1,0.4,8);
parm = 0.5/sqrt(k);
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'parm','M','N','method'};

N = 10;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
method = {'RT'};
% coreMethod = 'hp_FEM';
formulation = '';
M = 1;
plot3Dgeometry = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 6:8;
r = 10;

loopParameters = {'parm','M','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
M = 5:7;
plot3Dgeometry = 0;
M = 1;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
f = 1e3;
solveForPtot = true;
plot3Dgeometry = 1;
method = {'BEM'};
formulation = {'CCBIE','CBM'};
formulation = {'CBM'};
M = 5:7;
M = 1;
loopParameters = {'formulation','M','method','f'};
collectIntoTasks
