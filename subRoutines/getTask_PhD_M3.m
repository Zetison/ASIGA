scatteringCase = 'MS';

model = 'M3';  % Spherical shell

coreMethod = 'IGA';


f = 1e3;

M = 1:5;
parm = 1;
alpha = (0:0.5:180)*pi/180;
% alpha = (0:10:180)*pi/180;

plot2Dgeometry = 0;
plot3Dgeometry = 0;
degree = 2;
calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
applyLoad = 'planeWave';

loopParameters = {'M','N','formulation','method','f'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
M = 1:5;
method = {'IE'};
formulation = {'BGU'};
N = [1,3,5];
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
M = 1:6;
method = {'IENSG'};
N = [1,3,5,7,9];
formulation = {'BGU'};
collectIntoTasks