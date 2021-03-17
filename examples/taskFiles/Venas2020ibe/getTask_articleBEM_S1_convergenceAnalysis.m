misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'S1';

parm = 1;
misc.method = {'BEM'};
BC = 'SHBC';
% misc.coreMethod = {'IGA'};
misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
formulation = 'CCBIE';
% formulation = 'CRCBIE1';
c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
M = 1:7;
% M = 1:3;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
solveForPtot = true;

degree = 2;
calculateFarFieldPattern = 0;
err.calculateSurfaceError = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
solveForPtot = true;
loopParameters = {'M','misc.method','degree','misc.coreMethod','extraGPBEM','extraGP'};
% agpBEM = 2;
colBEM_C0 = 0;
% colBEM_C0 = 1/2;
quadMethodBEM = 'Adaptive';
% quadMethodBEM = 'Simpson';

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misc.coreMethod = {'IGA'};
degree = 3:4;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misc.coreMethod = {'linear_FEM'};
degree = 1;
M = [M, M(end)+1];

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formulation = 'SL2E';
misc.method = {'BA'};
misc.coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
degree = 2;
M = M(1:end-1);

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misc.coreMethod = {'IGA'};
degree = 3:4;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misc.coreMethod = {'linear_FEM'};
degree = 1;
M = [M, M(end)+1];

collectIntoTasks