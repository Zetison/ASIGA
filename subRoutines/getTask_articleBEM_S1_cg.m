scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';


% method = {'BA'};
% formulation = {'SL2E'};
method = {'BEM'};
formulation = {'GCBIE'};
coreMethod = {'IGA'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
% f = 1000;
% omega = 2*pi*f;
% k = omega/c_f;

M = 5;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

parm = 2;
degree = 4:5;
degree = 5;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

plotResidualError = true;
solveForPtot = true;
loopParameters = {'M','formulation','method','degree','coreMethod'};
collectIntoTasks
