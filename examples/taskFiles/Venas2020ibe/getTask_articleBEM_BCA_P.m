scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';
    
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

BC = 'NBC';

f = [1e2 1e3];             % Frequency
f = 1e2;             % Frequency
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = true;
calculateFarFieldPattern = true;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = 0;
applyLoad = 'pointPulsation';
method = {'BEM'};
formulation = {'CCBIE'};
M = [1,2];
% M = 1;
extraGPBEM = 50; % extra quadrature points around singularities for BEM formulations
storeSolution = 0;
storeFullVarCol = 0;
degree = [2,5];
% degree = 2;
agpBEM = 0.6;
% degree = 2;
solveForPtot = false;
loopParameters = {'method','degree','formulation','M','f'};
% quadMethodBEM = 'Adaptive2';
% collectIntoTasks

formulation = {'CBM','CCBIE'};
% formulation = {'CCBIEC'};
% degree = 5;
% M = 2;
degree = 2;
M = 1;
f = 1e3; 
% collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
f = 1e2;             % Frequency
collectIntoTasks




