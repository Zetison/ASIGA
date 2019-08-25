scatteringCase = 'BI';

model = 'S1';  % Spherical shell

coreMethod = {'IGA'};


f = 1e3;

alpha_s = 0;
beta_s = 0;   

M = 1:6;
parm = 2;

plot2Dgeometry = 0;
plot3Dgeometry = 0;
degree = 2;
% calculateVolumeError = 1;
calculateSurfaceError = 1;
computeCondNumber = false;
calculateFarFieldPattern = 0;
applyLoad = 'planeWave';
N = 1:2;

loopParameters = {'M','N','degree','formulation','coreMethod','method'};

%% IE simulation
method = {'IE'};
formulation = {'BGU'};
N = 1:3;
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = {'PS'};
M = 1:5;
delta = 0.5;
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEM simulation
method = {'IE'};
coreMethod = {'SEM'};
formulation = {'BGU'};
M = 1;
degree = 1:5;
N = 3;
collectIntoTasks


