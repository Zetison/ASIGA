scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S1';

coreMethod = 'IGA';
method = {'BEM'};
formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3','GCBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
formulation = {'CCBIE','CBM','GCBIE','GBM'};

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 

Mend = 5;
M = [4,Mend];
% M = 5;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;

parm = 1;
% degree = [4,5];
degree = [2,5];
runTasksInParallel = 0;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateSurfaceError = 1;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
solveForPtot = true;

loopParameters = {'extraGPBEM','agpBEM','M','degree','formulation','method','quadMethodBEM'};

agpBEM = 0.2:0.1:0.7;
extraGPBEM = 100; % extra quadrature points
quadMethodBEM = {'Adaptive'};
collectIntoTasks


extraGPBEM = [0,1,2,3,4,10:10:100]; % extra quadrature points
agpBEM = 0.7;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = {'BA'};
extraGPBEM = [extraGPBEM(1),extraGPBEM(end)];
agpBEM = [0.2,0.7];
formulation = {'SL2E'};
coreMethod = 'IGA';
collectIntoTasks