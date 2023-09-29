

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'S1';

f = 3e3;             % Frequency

parm = 2;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
beta = 30*pi/180;
calculateFarFieldPattern = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
misc.coreMethod = {'IGA'};
formulation = 'SL2E';
M = 1:2;
prePlot.plot3Dgeometry = 0;
degree = 2:7;
err.calculateSurfaceError = 1;
computeCondNumber = false;
solveForPtot = true;
loopParameters = {'M','misc.coreMethod','degree','misc.method'};
% collectIntoTasks

misc.applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries
misc.model = 'S1_P2';
% collectIntoTasks


plotResultsInParaview = 1;	% Only if misc.scatteringCase == 'Bi'
degree = 7;
M = 2;
collectIntoTasks

misc.applyLoad = 'planeWave'; % with analytic solution for arbitrary geometries
misc.model = 'S1';
collectIntoTasks
