scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

method = {'BEM'};
BC = 'SHBC';
coreMethod = 'IGA';
formulation = {'CCBIE','CBM'};
c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
M = 1:7;
% M = 6;
parm = [1,2];
parm = 1;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
solveForPtot = true;

colMethod = {'CG','Grev','GL'};
colMethod = {'Grev'};
degree = 2:5;
% degree = 3;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
solveForPtot = true;
loopParameters = {'colMethod','M','method','formulation','degree','parm'};
% agpBEM = 2;
colBEM_C0 = 0;
quadMethodBEM = 'Adaptive';

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formulation = {'SL2E'};
method = {'BA'};
colMethod = NaN;

collectIntoTasks