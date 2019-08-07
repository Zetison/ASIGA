getDefaultTaskValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'S1_P';  % Pulsating sphere
BC = 'NBC';

c_f = 1500;
k = 1;
omega = k*c_f;
f = omega/(2*pi); 

applyLoad = 'radialPulsation';

M = 4; 
parm = 1;

degree = 2;

plotFarField          = 0;
calculateSurfaceError = 0;
calculateFarFieldPattern = 0;
plot3Dgeometry = 0;
solveForPtot = true;
loopParameters = {'extraGPBEM','extraGP','agpBEM','colBEM_C0','method','formulation','quadMethodBEM'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
% formulation = {'CCBIE','CBM','GCBIE','GBM'};
formulation = {'CCBIE'};
quadMethodBEM = {'Simpson'};
colBEM_C0 = 0;
% colBEM_C0 = Inf;
extraGP = 0; % extra quadrature points
extraGPBEM = 8; % extra quadrature points around singularities for BEM formulations
% extraGPBEM = 12; % extra quadrature points around singularities for BEM formulations
% agpBEM = 2;
agpBEM = 2;

% collectIntoTasks

quadMethodBEM = {'Adaptive2'};
agpBEM = 1;
collectIntoTasks

