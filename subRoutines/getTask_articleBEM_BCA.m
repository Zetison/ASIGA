scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA'; % BeTSSi submarine
coreMethod = 'IGA';

% alpha_s = 270*pi/180;
alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

% f = [5e2, 1e3]; %[1e2 5e2 1e3];             % Frequency
f = [1e2 1e3];
f = 1e2;
alpha = (0:0.05:360)*pi/180;

plotResultsInParaview = 0;
plotMesh              = 0;	% Create additional Paraview files to visualize IGA mesh
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
BC = 'SHBC';

method = {'BEM'};
formulation = {'CCBIE'};
% formulation = {'CBM'};
M = 1:3;
degree = [2,5];
storeSolution = 0;
storeFullVarCol = 0;
agpBEM = 0.6;
loopParameters = {'method','formulation','M','degree','f','scatteringCase'};

collectIntoTasks

f = 1e3;
formulation = {'CBM'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calculateSurfaceError = 1;
applyLoad = 'radialPulsation';
BC = 'NBC';
M = 1:3;
degree = [2,5];
method = {'BA'};
formulation = {'SL2E'};
% collectIntoTasks

