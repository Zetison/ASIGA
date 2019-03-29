

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

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
method = {'BA'};
coreMethod = {'IGA'};
formulation = 'SL2E';
M = 1:2;
plot3Dgeometry = 0;
degree = 2:7;
calculateSurfaceError = 1;
computeCondNumber = false;
loopParameters = {'M','coreMethod','degree','method'};
% collectIntoTasks

applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
model = 'S1_P2';
% collectIntoTasks


plotResultsInParaview = 1;	% Only if scatteringCase == 'Bi'
degree = 7;
M = 2;
collectIntoTasks

applyLoad = 'planeWave'; % with analytic solution for arbitrary geometries
model = 'S1';
collectIntoTasks
