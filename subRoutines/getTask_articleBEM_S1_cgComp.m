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
M = 1:6;
parm = [0,1];

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
solveForPtot = true;

colMethod = {'CG','Grev','GL'};
degree = 2:5;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
loopParameters = {'colMethod','M','parm','method','formulation','degree'};
% agpBEM = 2;
colBEM_C0 = 0;
quadMethodBEM = 'Adaptive';

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formulation = {'SL2E'};
method = {'BA'};
colMethod = NaN;

collectIntoTasks