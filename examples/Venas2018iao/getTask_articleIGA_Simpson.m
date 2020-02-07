%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IE simulation
scatteringCase = 'BI';
model = 'SS';  % Spherical shell
% method = {'IE','IENSG'};
method = {'IE'};
BC = 'SHBC';
formulation = {'PGU','PGC','BGU','BGC'};
% formulation = {'BGU'};
coreMethod = 'IGA';
computeCondNumber = 0;
runTasksInParallel = 0;

c_f = 1500; % Speed of sound in outer fluid
% k = 100;             % Wave number for outer fluid domain
k = 2;             % Wave number for Simpson2014aib
omega = c_f*k;   % Angular frequency
f = omega/(2*pi);    % Frequency

M = 3;

alpha_s = pi;
beta_s = 0;  
alpha = (0:0.5:360)*pi/180;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometry  

plotFarField = false; 
r = 5; % radii for near-field evaluation
degree = 3;

N = 4;

loopParameters = {'M','method','formulation'};
calculateSurfaceError = 1;
parm = 1;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
useNeumanProj = 0;
solveForPtot = true;
formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
method = {'BEM'};
% formulation = {'CCBIE', 'CHBIE', 'CBM','GCBIE', 'GHBIE', 'GBM'};
formulation = {'CCBIE'};
% formulation = 'CBM';
% useNeumanProj = [0,1];
colBEM_C0 = 0;
solveForPtot = true;
% extraGP = 2; % extra quadrature points
% extraGPBEM = 14; % extra quadrature points around singularities for BEM formulations
collectIntoTasks
