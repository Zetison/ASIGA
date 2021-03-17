misc.scatteringCase = 'MS';

misc.model = 'M1';  % Spherical shell

misc.coreMethod = 'IGA';


f = [1e2,1e3];
f = 1e3;

parm = [];
alpha = (0:0.1:360)*pi/180;

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
degree = 2;
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
misc.applyLoad = 'planeWave';

loopParameters = {'M','N','formulation','misc.method','f'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
M = 4:6;
M = 4;
misc.method = {'IE'};
formulation = {'BGU'};
N = [3,5,7];
misc.method = {'IENSG'};
N = [3,5];
N = 5;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
M = 1;
N = NaN;
formulation = {'GBM'};
% formulation = {'CCBIE'};
solveForPtot = true;
% collectIntoTasks