scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

coreMethod = 'IGA';
method = {'BEM'};
formulation = {'CCBIE','CBM','GCBIE','GBM'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 

Mend = 4;
M = [1,Mend];
% M = Mend;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

parm = 1;
degree = [2,5];
runTasksInParallel = 1;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

extraGP = 0; % extra quadrature points
extraGPBEM = [0,1,2,4,8,16,32,64]; % extra quadrature points
agpBEM = 1:10;
colBEM_C0 = 1/2;
quadMethodBEM = {'New','Simpson'};

loopParameters = {'extraGPBEM','agpBEM','M','degree','formulation','method','quadMethodBEM'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'BA'};
formulation = {'SL2E'};
coreMethod = {'IGA'};
collectIntoTasks