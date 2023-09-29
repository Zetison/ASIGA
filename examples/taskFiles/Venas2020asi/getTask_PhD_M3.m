misc.scatteringCase = 'MS';

misc.model = 'M3';  % Spherical shell

misc.coreMethod = 'IGA';


f = 1e3;

M = 1:5;
parm = 1;
alpha = (0:0.1:180)*pi/180;
% alpha = (0:10:180)*pi/180;

prePlot.plot2Dgeometry = 1;
prePlot.plot3Dgeometry = 1;
degree = 2;
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
misc.applyLoad = 'planeWave';

loopParameters = {'M','N','formulation','misc.method','f'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
M = 1:5;
M = 3;
misc.method = {'IE'};
formulation = {'BGU'};
N = [1,3];
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
M = 1:6;
misc.method = {'IENSG'};
N = [1,3,5,7,9];
formulation = {'BGU','BGC'};
collectIntoTasks