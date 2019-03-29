

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';


method = {'BEM'};
formulation = {'GBM'};
formulation = {'CCBIE'};
coreMethod = {'SEM'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 
% f = 1000;
% omega = 2*pi*f;
% k = omega/c_f;

Mend = 3;
M = 1:Mend;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha = (0:0.5:360)*pi/180;
beta = 30*pi/180;

parm = 1;
degree = 2:5;
% degree = 2;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

extraGP = 0:4; % extra quadrature points
extraGPBEM = 0:6; % extra quadrature points
extraGP = -1;
extraGP = 0;
extraGPBEM = 4;
loopParameters = {'M','formulation','method','degree','coreMethod','extraGP','extraGPBEM'};
collectIntoTasks
% 
% coreMethod = {'linear_FEM'};
% degree = 1;
% M = 1:Mend+1;
% % collectIntoTasks
% 
% coreMethod = {'IGA'};
% M = 1:Mend;
% degree = 3:4;
% % collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = {'BA'};
formulation = {'SL2E'};
coreMethod = {'IGA'};
collectIntoTasks
% 
% coreMethod = {'linear_FEM'};
% degree = 1;
% M = 1:Mend+1;
% % collectIntoTasks
% 
% % degree = 3:4;
% degree = 2;
% % coreMethod = {'XI','IGA'};
% coreMethod = {'IGA'};
% M = 1:Mend;
% % M = 2;
% collectIntoTasks
% 
% 
% 
