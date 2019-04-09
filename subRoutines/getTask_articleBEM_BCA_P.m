scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'BCA_P'; % BeTSSi submarine
coreMethod = 'IGA';
    
plot3Dgeometry = 0;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometr

BC = 'SHBC';

f = 1e2; %[1e2 5e2 1e3];             % Frequency
plotResultsInParaview = 1;
plotMesh              = 1;	% Create additional Paraview files to visualize IGA mesh
calculateSurfaceError = true;
calculateFarFieldPattern = true;
plotTimeOscillation   = 0;	% Create 30 paraview files in order to visualize a dynamic result
computeCondNumber = true;
applyLoad = 'radialPulsation';
method = {'BEM'};
formulation = {'CCBIE'};
M = 1:2;
storeSolution = 0;
storeFullVarCol = 0;
degree = [2,6];
loopParameters = {'method','formulation','M','degree'};
collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
collectIntoTasks




