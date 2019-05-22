
% scatteringCase = 'MS';  % Monostatic scattering
scatteringCase = 'BI'; % Bistatic scattering
% scatteringCase = 'Sweep'; % (monostatic) Frequency sweep

% model = 'PS';  % Pulsating sphere
model = 'SS';  % Spherical shell
% model = 'M1';  % BeTSSi model 1
% model = 'M3';  % BeTSSi model 3
% model = 'M4';  % BeTSSi model 4
% model = 'M5A'; % BeTSSi model 5A
% model = 'M5B'; % BeTSSi model 5B
% model = 'MS';  % Mock Shell
% model = 'BA';  % Barrel
% model = 'EL';  % Ellipsoid

method = 'IE';
% method = 'IENSG';
% method = 'BA'; % only when analytic solution exist
% method = 'BEM';

% Alternatives for infinite elements (IE and IENSG):
% formulation = 'PGU'; % - Petrov-Galerkin Unconjugated (unconjugated Lies IEM)
% formulation = 'PGC'; % - Petrov-Galerkin (conjugated Lies IEM)
% formulation = 'BGU'; % - Bubnov-Galerkin Unconjugated (unconjugated Burnett IEM)
% formulation = 'BGC'; % - Bubnov-Galerkin Conjugated (conjugated Burnett IEM)

% Alternatives for boundary element method:
% formulation = 'CCBIE'; % - Collocation Conventional Boundary Integral Equation
% formulation = 'CHBIE'; % - Collocation Hypersingular Boundary Integral Equation
% formulation = 'CBM';   % - Collocation Burton Miller
% formulation = 'GCBIE'; % - Galerkin Conventional Boundary Integral Equation
% formulation = 'GHBIE'; % - Galerkin Hypersingular Boundary Integral Equation
% formulation = 'GBM';   % - Galerkin Burton Miller

% Alternatives for best approximation (BA):
% formulation = 'Surf';



BC = 'SHBC';
plot2Dgeometry        = false;  % Plot cross section of mesh and geometry
plot3Dgeometry        = false;  % Plot visualization of mesh and geometry in 3D
calculateSurfaceError = false;	% Only if scatteringCase == 'Bi'
calculateVolumeError  = false;	% Only if scatteringCase == 'Bi'
plotResultsInParaview = false;	% Only if scatteringCase == 'Bi'
plotTimeOscillation   = false;	% Create 30 paraview files in order to visualize a dynamic result
plotMesh              = false;	% Create additional Paraview files to visualize IGA mesh
plotFarField          = true;	% If false, plots the near-field instead
useNeumanProj         = false;
runTasksInParallel    = false;
calculateFarFieldPattern = true;
computeCondNumber = false;
storeSolution = false;
storeFullVarCol = false;
clearGlobalMatrices = true;
extraGP = 0; % extra quadrature points
extraGPBEM = 4; % extra quadrature points around singularities for BEM formulations
agpBEM = 2; % parameter for adaptiv Gauss point integration around singularities for BEM formulations
useROM = false;
exteriorProblem = true;
initMeshFactXi = 1;
initMeshFactZeta = 1;
colBEM_C0 = 0;
quadMethodBEM = 'New';
M = 1;

applyLoad = 'planeWave';
% applyLoad = 'radialPulsation'; % with analytic solution for arbitrary geometries
% applyLoad = 'SimpsonTorus'; % with analytic solution for arbitrary geometries

coreMethod = 'IGA';

alpha_s = NaN;
beta_s = NaN;

alpha = (0:0.5:360)*pi/180;
beta = 0;

r = 1; % radii for near-field evaluation. Assume by default plotFarField = true

N_max = inf; % number of terms in analytic solution for scattering on spherical shell


loopParameters = {'M'};


