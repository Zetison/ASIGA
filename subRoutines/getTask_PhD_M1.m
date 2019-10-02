scatteringCase = 'MS';

model = 'M1';  % Spherical shell

coreMethod = 'IGA';


f = [1e2,1e3];
f = 1e3;

parm = [];
alpha = (0:0.1:360)*pi/180;

plot2Dgeometry = 0;
plot3Dgeometry = 0;
degree = 2;
calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
applyLoad = 'planeWave';

loopParameters = {'M','N','formulation','method','f'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
M = 4:6;
M = 4;
method = {'IE'};
formulation = {'BGU'};
N = [3,5,7];
method = {'IENSG'};
N = [3,5];
N = 5;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
M = 1;
N = NaN;
formulation = {'GBM'};
% formulation = {'CCBIE'};
solveForPtot = true;
% collectIntoTasks