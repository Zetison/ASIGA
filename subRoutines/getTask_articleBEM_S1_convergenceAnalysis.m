scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

parm = 1;
method = {'BEM'};
BC = 'SHBC';
% coreMethod = {'IGA'};
coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
formulation = 'CCBIE';
c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
M = 1:7;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

degree = 2;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
loopParameters = {'M','method','degree','coreMethod','extraGPBEM','extraGP'};
extraGPBEM = 0;
extraGP = 0;

% collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'IGA'};
degree = 3:4;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'linear_FEM'};
degree = 1;
M = [M, M(end)+1];

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formulation = 'SL2E';
method = {'BA'};
coreMethod = {'hp_FEM','h_FEM','C0_IGA','IGA'};
degree = 2;
M = 1:7;

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'IGA'};
degree = 3:4;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'linear_FEM'};
degree = 1;
M = [M, M(end)+1];

collectIntoTasks