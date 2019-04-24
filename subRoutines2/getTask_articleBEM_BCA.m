

scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
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
% BC = 'SSBC';
BC = 'SHBC';

method = {'BEM'};
% coreMethod = 'XI';
formulation = {'CCBIE'};
% formulation = {'CRCBIE'};
M = 3;
degree = 5;
storeSolution = 0;
storeFullVarCol = 0;
loopParameters = {'method','formulation','M','degree','f','formulation','scatteringCase'};
extraGPBEM = 0;
extraGP = 0;
collectIntoTasks

M = 1:2;
degree = [2,4,6];
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 1e3; %[1e2 5e2 1e3];             % Frequency
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = 0;
calculateFarFieldPattern = 0;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = 0;
applyLoad = 'radialPulsation';
model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';
method = {'BEM'};
formulation = 'CCBIE';
method = {'BA'};
formulation = 'SL2E';
M = 1;
storeSolution = 0;
storeFullVarCol = 0;
degree = 5;
% degreeElev = 0:8;
loopParameters = {'method','degree'};
% collectIntoTasks




