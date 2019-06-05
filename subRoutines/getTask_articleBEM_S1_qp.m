scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

coreMethod = 'IGA';
method = {'BEM'};
formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3','GCBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
formulation = {'GCBIE','GBM'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 

Mend = 5;
M = [1,Mend];
M = 5;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

parm = 1;
% degree = [4,5];
degree = 4;
runTasksInParallel = 0;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

extraGP = 0; % extra quadrature points
extraGPBEM = [0,4,8,16,32,64]; % extra quadrature points
agpBEM = 0:0.5:1.5;
% agpBEM = 3;
colBEM_C0 = 0;
quadMethodBEM = {'Simpson'};

loopParameters = {'extraGPBEM','agpBEM','M','degree','formulation','method','quadMethodBEM'};
% collectIntoTasks

agpBEM = 0.2:0.4:1.4;
quadMethodBEM = {'Adaptive'};
collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'BA'};
extraGPBEM = NaN;
agpBEM = [agpBEM(1),agpBEM(end)];
formulation = {'SL2E'};
coreMethod = 'IGA';
collectIntoTasks