scatteringCase = 'BI';      % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
model = 'S1';               % Geometry model
f = 1e3; % Frequency
M = 1:3;                    % Mesh
alpha_s = 0;                % Aspect angle for source
beta_s = 0;                 % Elevation angle for source
alpha = (0:0.5:360)*pi/180; % Aspect angle for receiver
beta = beta_s;              % Elevation angle for receiver
calculateVolumeError = 0;   
calculateSurfaceError = 1;
calculateFarFieldPattern = 1;
plot3Dgeometry = 0;
degree = 4;
parm = 2;                   % use 2nd parametrization of the sphere (6 patches)
coreMethod = 'IGA';         % Underlying mesh type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
method = {'IE'};
formulation = 'BGU';
N = 6;
IEbasis = 'Chebyshev';
loopParameters = {'M','parm','method'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = 'PS'; % PS = point solution, SS = spherical solution
computeCondNumber = false;
delta = 0.2;
extraGP = 2;
loopParameters = {'M','delta','extraGP','method'};

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = 'SL2E';
plot3Dgeometry = 0;
loopParameters = {'M','parm','method'};
collectIntoTasks

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
formulation = 'CCBIE';
plot3Dgeometry = 0;
degree = 4;
loopParameters = {'M','method'};
solveForPtot = true;
collectIntoTasks

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
formulation = '';
plot3Dgeometry = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
collectIntoTasks


