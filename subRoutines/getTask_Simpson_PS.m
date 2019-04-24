getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'S1_P';  % Pulsating sphere
method = {'IE','IENSG'};
BC = 'NBC';
formulation = {'BGU'};
coreMethod = 'IGA';

f = 1/(2*pi);

applyLoad = 'radialPulsation';

M = 1:3; 
M = 1; 
parm = 1;

degree = 4;
N = 6;

plotFarField          = 0;
calculateSurfaceError = 1;
calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
loopParameters = {'M', 'method','formulation'};

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = {'SL2E'};
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
method = {'BEM'};
formulation = {'CCBIE'};
collectIntoTasks
