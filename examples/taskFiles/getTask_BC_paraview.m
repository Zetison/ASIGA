

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BC'; % BeTSSi submarine
coreMethod = 'IGA';
% coreMethod = {'IGA'};

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

f = 1e3;
method = 'IE';
formulation = 'BGU';
IEbasis = 'Lagrange';
M = 1;
degreeElev = 0;
plotResultsInParaview = 1;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
% BC = 'SSBC';
BC = 'SHBC';
N = 3;
% loopParameters = {'M','degreeElev','f'};
loopParameters = {'f','BC'};

collectIntoTasks