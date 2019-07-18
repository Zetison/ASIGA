scatteringCase = 'BI';

model = 'S1';  % Spherical shell

coreMethod = 'IGA';


f = 20e3;

alpha_s = 0;
beta_s = 0;   

M = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
method = {'RT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degree = 2;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 1;
N = 3:6;

loopParameters = {'N','M','method'};
collectIntoTasks
