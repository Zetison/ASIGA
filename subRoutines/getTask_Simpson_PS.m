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
%% BEM simulation
runTasksInParallel = 1;
method = {'BEM'};
% formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE','CHBIE','CBM','CRCBIE1','CRCBIE2','CRCBIE3'};
% formulation = {'CCBIE'};
% formulation = {'CCBIE','CBM'};
% formulation = {'CBM'};
% formulation = {'CHBIE'};

colBEM_C0 = [1/2,0]; 
% colBEM_C0 = 0;
% extraGP = [0,1,2,4,8,16]; % extra quadrature points
extraGP = 0; % extra quadrature points
extraGPBEM = 50; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = [30,40,50,60]; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 10; % extra quadrature points around singularities for BEM formulations
% agpBEM = [2,4,8,16];
% agpBEM = 10;
% useNeumanProj = [0,1];
agpBEM = 1:12;
quadMethodBEM = {'Simpson'};

loopParameters = {'formulation','extraGPBEM','extraGP','agpBEM','colBEM_C0','method','quadMethodBEM'};

collectIntoTasks

agpBEM = (1:12)/5;
% agpBEM = 2;
% colBEM_C0 = 1/2; 
quadMethodBEM = {'Adaptive'};
% quadMethodBEM = {'Simpson'};
collectIntoTasks

agpBEM = 1:12;
formulation = {'GCBIE','GHBIE','GBM','GRCBIE1','GRCBIE2','GRCBIE3'};
quadMethodBEM = {'Simpson'};
colBEM_C0 = NaN;
collectIntoTasks


agpBEM = (1:12)/5;
quadMethodBEM = {'Adaptive'};
colBEM_C0 = NaN;
collectIntoTasks

