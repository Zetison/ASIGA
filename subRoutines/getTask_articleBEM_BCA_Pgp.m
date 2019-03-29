

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = 1e2; %[1e2 5e2 1e3];             % Frequency
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
calculateSurfaceError = true;
calculateFarFieldPattern = false;

BC = 'SHBC';

applyLoad = 'radialPulsation';
method = {'BEM'};
formulation = 'CCBIE';
M = 1;
storeSolution = false;
storeFullVarCol = false;
loopParameters = {'method','degree','extraGP','extraGPBEM','agpBEM'};
degree = 4;

extraGP = 0; % extra quadrature points
extraGPBEM = 0; % extra quadrature points
agpBEM = [0.25,0.5,1,1.5,2,2.5,3];
collectIntoTasks

% extraGP = [0,4];
% extraGPBEM = 0;
% method = {'BA'};
% formulation = 'SL2E';
% collectIntoTasks
