getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'SS_P';  % Pulsating sphere
method = {'IENSG','IE'};
BC = 'SHBC';
formulation = 'PGU';
coreMethod = 'IGA';

f = 1/(2*pi);

applyLoad = 'radialPulsation';

M = 1:3; 
M = 1; 
parm = 2;

degreeElev = 1;
N = 1;

plotFarField          = 1;
calculateSurfaceError = 1;
calculateVolumeError  = 1;
plot3Dgeometry = 0;
loopParameters = {'M', 'method'};

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = 'SL2E';
calculateVolumeError  = 0;
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RBEM simulation
method = {'BEM'};
formulation = 'CRCBIE';
M = 1; 
plotMesh = 0;
plotResultsInParaview = 0;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
formulation = 'CCBIE';
collectIntoTasks
