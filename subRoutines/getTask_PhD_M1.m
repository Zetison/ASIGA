scatteringCase = 'MS';

model = 'M1';  % Spherical shell

coreMethod = 'IGA';


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
M = 4:5;
% M = 2;
method = {'IE'};
formulation = {'BGU'};
N = [3,5,7];
method = {'IENSG'};
% N = [1,3,5,7,9];
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
% M = 2;
N = NaN;
formulation = {'CBM','GBM'};
% formulation = {'CCBIE'};
solveForPtot = true;
collectIntoTasks