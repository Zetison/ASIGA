misc.scatteringCase = 'BI';

misc.model = 'S1';  % Spherical shell

misc.coreMethod = 'IGA';


f = 1e3;

alpha_s = 0;
beta_s = 0;   

M = 1:6;
parm = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC simulation
misc.method = {'ABC'};
% formulation = {'HH','BGT'};
formulation = {'BGT'};
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
degree = 2;
% calculateVolumeError = 1;
err.err.calculateSurfaceError = 1;
computeCondNumber = false;
calculateFarFieldPattern = 0;
misc.applyLoad = 'planeWave';
N = 1:2;

loopParameters = {'M','N','degree','formulation','misc.method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
misc.method = {'IE'};
formulation = {'BGC','BGU','PGU','PGC'};
N = [1,3,6,9];
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IENSG simulation
misc.method = {'IENSG'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
formulation = {'SL2E'};
N = NaN;
collectIntoTasks


