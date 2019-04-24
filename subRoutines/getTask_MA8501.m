

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = 'IE';
formulation = 'PGU';
% coreMethod = {'SEM','IGA'};
coreMethod = {'SEM'};
% coreMethod = {'IGA'};
% f = [1e3,1e4];             % Frequency
f = 1e3;             % Frequency
% f = 3.98e4;             % Frequency


M = 1; %3:6
N = [3,9];
N = 3;
N_max = N-1;
% alpha_s = 0;
% beta_s = pi/2;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = pi;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;
calculateVolumeError = 1;
calculateSurfaceError = 0;
computeCondNumber = true;
clearGlobalMatrices = false;
% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries

calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
degree = 0:15;
degree = 0:5;
parm = 2;

loopParameters = {'degree','N','coreMethod','f'};
collectIntoTasks
