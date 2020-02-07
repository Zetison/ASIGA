

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
BC = 'NBC';

method = {'BEM'};
coreMethod = 'IGA';
formulation = 'CCBIE';
M = 1;
degree = 0:10;
loopParameters = {'method','degreeElev'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = true;
calculateFarFieldPattern = false;
computeCondNumber = false;
applyLoad = 'pointPulsation';
model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';
method = {'BA','BEM'};
formulation = 'SL2E';
M = 1;
degree = 0;
% degreeElev = 0:8;
loopParameters = {'method','degreeElev'};
collectIntoTasks
