getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'S1_P';  % Pulsating sphere
method = {'IE'};
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
initMeshFactZeta = 32;

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
runTasksInParallel = 1;
method = {'BEM'};
% formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE','CBM','CRCBIE3'};
% formulation = {'CCBIE'};
formulation = {'CCBIE','CBM'};
% formulation = {'CBM'};
% formulation = {'CHBIE'};

colBEM_C0 = [1/2,0];
% colBEM_C0 = 2;
% extraGP = [0,1,2,4,8,16]; % extra quadrature points
extraGP = 0; % extra quadrature points
extraGPBEM = 50; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 0; % extra quadrature points around singularities for BEM formulations
% agpBEM = [2,4,8,16];
agpBEM = 1:10;
% agpBEM = 2;
quadMethodBEM = {'Simpson'};

loopParameters = {'extraGPBEM','extraGP','agpBEM','colBEM_C0','method','formulation','quadMethodBEM'};

collectIntoTasks

% agpBEM = [0.6,1.2,2.4];
agpBEM = 1:6;
quadMethodBEM = {'New'};
collectIntoTasks

formulation = {'GCBIE','GBM'};
agpBEM = 1:10;
quadMethodBEM = {'Simpson'};
colBEM_C0 = NaN;
collectIntoTasks

formulation = {'GCBIE','GBM'};
agpBEM = 1:6;
quadMethodBEM = {'New'};
colBEM_C0 = NaN;
collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
formulation = {'CRCBIE1','CRCBIE2'};
extraGP = 0; % extra quadrature points
extraGPBEM = 0:32; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 32; % extra quadrature points around singularities for BEM formulations
agpBEM = 0;
colBEM_C0 = 2;
% collectIntoTasks

