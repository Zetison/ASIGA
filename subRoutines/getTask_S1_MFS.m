% Figure 6 in Michler2007itp (page 842)

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = 'MFS';
% formulation = {'BGU','PGU','BGC','PGC'};
formulation = '';
% formulation = {'BGU'};
% coreMethod = {'C0_IGA'};
coreMethod = 'IGA';
% coreMethod = {'hp_FEM'};
k = 1;
c = 1500;
omega = c*k;
f = omega/(2*pi);

M = 1; %3:6
% N = [3,9];
% N = 3;
alpha_s = 0;
beta_s = pi/2;

calculateVolumeError = 0;
calculateSurfaceError = 1;

calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
plot2Dgeometry = 0;
degree = 4;
% parm = linspace(0.25,0.35,10);
parm = 0.32778;
loopParameters = {'M','parm'};
collectIntoTasks

