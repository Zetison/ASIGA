misc.scatteringCase = 'BI';

misc.model = 'Shirron';  % Spherical shell

misc.coreMethod = 'IGA';


c_f = 1500;  % Speed of sound in fluid domains
k = 10;
omega = k*c_f;
f = omega/(2*pi);

parm = [];
alpha = (0:0.1:180)*pi/180;
% alpha = (0:10:180)*pi/180;
alpha_s = [180,0]*pi/180;
alpha_s = 0;
beta_s = 0;
parm = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 1;
degree = 3;
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
calculateFarFieldPattern = 1;
misc.applyLoad = 'planeWave';

loopParameters = {'M','N','formulation','misc.method','f','alpha_s'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
M = 1;
misc.method = {'IE'};
formulation = {'BGU'};
N = [3,5,7];
N = 3;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
% M = 5;
% N = 3;
misc.method = {'IENSG'};
N = 3;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
% M = 1;
% M = 5;
N = NaN;
formulation = {'GBM'};
solveForPtot = true;
collectIntoTasks