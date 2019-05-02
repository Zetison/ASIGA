getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'S1_P';  % Pulsating sphere
method = {'IE','IENSG'};
BC = 'NBC';
formulation = {'BGU'};
coreMethod = 'IGA';


c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 

applyLoad = 'radialPulsation';

M = 1; 
parm = 2;

degree = 4;
N = 6;
extraGP = 0:2; % extra quadrature points

plotFarField          = 0;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
loopParameters = {'extraGPBEM','extraGP','agpBEM','colBEM_C0','method','formulation'};

% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RBEM simulation
method = {'BEM'};
formulation = {'CRCBIE1','CRCBIE2','CRCBIE3'};
extraGP = 2; % extra quadrature points
extraGPBEM = 2; % extra quadrature points around singularities for BEM formulations
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
runTasksInParallel = 0;
method = {'BEM'};
% formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE','CBM'};
formulation = {'CRCBIE1','CRCBIE2','CRCBIE3'};
colBEM_C0 = [2,Inf];
% colBEM_C0 = Inf;
% extraGP = [0,1,2,4,8,16]; % extra quadrature points
extraGP = [0,16]; % extra quadrature points
extraGPBEM = 0:32; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 32; % extra quadrature points around singularities for BEM formulations
% agpBEM = [2,4,8];
agpBEM = [2,8];
collectIntoTasks

formulation = {'GCBIE','GBM'};
colBEM_C0 = Inf;
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = {'SL2E'};
extraGP = 16; % extra quadrature points
extraGPBEM = [extraGPBEM(1),extraGPBEM(end)]; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 32; % extra quadrature points around singularities for BEM formulations
agpBEM = NaN;
% collectIntoTasks

