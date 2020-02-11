scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'S5';

parm = 2;
plot2Dgeometry = 0;
plot3Dgeometry = 1;
method = {'BEM'};
formulation = 'GCBIE';
solveForPtot = true;
% BC = 'NNBC';
% BC = 'SSBC';
BC = 'SHBC';
f = 1e1; 
M = 1:2;
% M = 1:3; %1:5
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha_s = 0;
beta_s = pi/2;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plotResultsInParaview = 1;
initMeshFactZeta = 1;
applyLoad = 'radialPulsation';
% applyLoad = 'planeWave';

extraGP = 0; % extra quadrature points
degree = 4;
if parm == 2 && degree < 4
    degree = 4;
end
degree = 2;
loopParameters = {'M','method'};
collectIntoTasks

method = {'BA'};
formulation = 'SL2E';
solveForPtot = false;
% formulation = 'VL2E';
collectIntoTasks
    
method = {'IE'};
formulation = 'BGU';
solveForPtot = false;
calculateVolumeError = 0;
collectIntoTasks
    