

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';

alpha_s = 270*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = 2e2; %[1e2 5e2 1e3];             % Frequency
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
% BC = 'SSBC';
BC = 'SHBC';

method = {'BEM'};
coreMethod = 'IGA';
formulation = 'CCBIE';
M = 1;
degree = 0:10;
loopParameters = {'method','degreeElev'};
collectIntoTasks