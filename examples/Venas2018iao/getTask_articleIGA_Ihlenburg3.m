

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';

method = {'IE'};
formulation = 'BGU';

BC = 'NNBC';

coreMethod = 'IGA';

c_f = 1524;
k = 2;
omega = k*c_f;
f = omega/(2*pi); 
parm = 1;


M = 5;

N = 6;

alpha_s = pi;
beta_s = 0;

degree = 2;
plotResultsInParaview = 1;
calculateFarFieldPattern = 0;
calculateVolumeError = 1;
calculateSurfaceError = 0;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
scaleForErrorPlot = 0;
plotMesh = 1;
loopParameters = {'M','method'};

collectIntoTasks
    
method = {'BA'};
formulation = 'VL2E';
collectIntoTasks
