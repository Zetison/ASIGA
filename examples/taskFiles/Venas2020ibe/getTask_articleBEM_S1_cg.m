misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'S1';


% misc.method = {'BA'};
% formulation = {'SL2E'};
misc.method = {'BEM'};
formulation = {'GCBIE'};
misc.coreMethod = {'IGA'};

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
err.calculateSurfaceError = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;

plotResidualError = true;
solveForPtot = true;
loopParameters = {'M','formulation','misc.method','degree','misc.coreMethod'};
collectIntoTasks
