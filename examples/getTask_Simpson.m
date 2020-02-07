%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'SS';  % Spherical shell
% method = {'IE','IENSG'};
method = {'IE'};
BC = 'SHBC';
% formulation = {'PGU','PGC','BGU','BGC'};
formulation = {'BGU'};
coreMethod = 'IGA';
computeCondNumber = 0;
runTasksInParallel = 0;

c_f = 1500; % Speed of sound in outer fluid
% k = 100;             % Wave number for outer fluid domain
k = 2;             % Wave number for Simpson2014aib
omega = c_f*k;   % Angular frequency
f = omega/(2*pi);    % Frequency

M = 1:3;

alpha_s = pi;
beta_s = 0;  
alpha = (0:0.1:360)*pi/180;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometry  

if 0
    plotFarField = true;
    calculateSurfaceError = true;
    LpOrder = 2; % For error calculation in calcSurfError()
    calculateVolumeError  = true;
    degree = 2;
else % reproduce plot in Simpson2014aib
    plotFarField = false; 
    r = 5; % radii for near-field evaluation
    degree = 3;
end
N = 4;
useNeumanProj = 0;

loopParameters = {'M','method','formulation','useNeumanProj'};
parm = 1;
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
useNeumanProj = 0;
formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ABC simulation
method = {'ABC'};
formulation = 'HH'; % Hagstrøm Harian
N = 1:2;
loopParameters = {'M','method','N'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
M = 1:6;
degree = 2;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
% collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
% formulation = {'CCBIE', 'CHBIE', 'CBM'};
formulation = {'CCBIE'};
% formulation = 'CBM';
useNeumanProj = [0,1];
% colBEM_C0 = 2;
% extraGP = 2; % extra quadrature points
% extraGPBEM = 14; % extra quadrature points around singularities for BEM formulations
solveForPtot = true;
collectIntoTasks