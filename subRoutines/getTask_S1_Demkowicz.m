% Figure 6 in Michler2007itp (page 842)

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'IE'};
% formulation = {'BGU','PGU','BGC','PGC'};
formulation = {'BGU','BGC'};
% formulation = {'BGU'};
% coreMethod = {'C0_IGA'};
coreMethod = {'IGA','hp_FEM','C0_IGA','SEM'};
% coreMethod = {'hp_FEM'};
lambda = 1;
r_a = 2; % radius of artificial boundary
k = 2*pi/lambda;
c = 1500;
omega = c*k;
f = omega/(2*pi);
% f = 1e3;             % Frequency
% f = 3.98e4;             % Frequency

M = 2; %3:6
% N = [3,9];
% N = 3;
alpha_s = 0;
beta_s = pi/2;

calculateVolumeError = 1;
calculateSurfaceError = 0;

calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
plot2Dgeometry = 0;
degree = 4:9;
% degree = 7;
parm = [1,2];
parm = 1;

loopParameters = {'degree','M','formulation','coreMethod','method'};
collectIntoTasks

method = {'BA'};
formulation = {'VL2E'};
coreMethod = {'hp_FEM','C0_IGA'};
collectIntoTasks

method = {'ABC'};
coreMethod = {'C0_IGA'};
N = [1,2];
loopParameters = {'degree','M','N','formulation','coreMethod','method'};
formulation = {'BGT'};
collectIntoTasks


