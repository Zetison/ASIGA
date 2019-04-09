

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BC_P'; % BeTSSi submarine
coreMethod = 'IGA';
% coreMethod = {'IGA'};

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

f = 5e2;
applyLoad = 'radialPulsation';
method = 'IE';
formulation = 'BGU';
IEbasis = 'Lagrange';
M = 2;
degreeElev = 1;
calculateSurfaceError = 1;
calculateVolumeError  = 1;
plotResultsInParaview = 1;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
BC = 'NBC';
N = 3;

collectIntoTasks