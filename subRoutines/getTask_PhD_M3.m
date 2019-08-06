scatteringCase = 'MS';

model = 'M3';  % Spherical shell

coreMethod = 'IGA';


f = 1e2;

M = 1:3;
parm = 1;
alpha = (0:0.5:180)*pi/180;
% alpha = (0:10:180)*pi/180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC simulation
method = {'ABC'};
% formulation = {'HH','BGT'};
formulation = {'BGT'};
plot2Dgeometry = 0;
plot3Dgeometry = 0;
degree = 2;
% calculateVolumeError = 1;
calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
applyLoad = 'planeWave';
N = 1:2;

loopParameters = {'M','N','formulation','method'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
method = {'IE'};
formulation = {'BGU'};
N = 1:3;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
method = {'IENSG'};
formulation = {'BGU'};
collectIntoTasks