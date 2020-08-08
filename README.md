# ASIGA
The ASIGA (Acoustic Scattering using IsoGeometric Analysis) toolbox provides a framework for simulating acoustic scattering problems using IGA. 

## Installation
Get the toolbox from GitHub (the e3Dss repository is needed for exact solutions and export_fig are needed for beautiful graphics for publications):
```
git clone https://github.com/altmany/export_fig
git clone https://github.com/Zetison/e3Dss 
git clone https://github.com/Zetison/ASIGA
```

## Run program
The following will run the input script scriptName.m
```
matlab -nodisplay -nodesktop -nosplash -r "cd ASIGA, studyName = 'scriptName'; main; quit"
```
Alternatively, from the MATLAB GUI a file availableStudies.m must be created in the parent folder of the ASIGA folder which should contain the command 

```Matlab
studyName = '<scriptName>';
```

## Models
The following models are available
- Sphere or spherical shells (i.e. model = 'S1')
- BeTSSi model 1 (model = 'M1')
- BeTSSi model 2 (model = 'M2')
- BeTSSi model 3 (model = 'M3')
- BeTSSi model 4 (model = 'M4')
- BeTSSi model 5A (model = 'M5A')
- BeTSSi model 5B (model = 'M5B')
- BeTSSi pressure hull (model = 'PH')
- Stripped BeTSSi submarine (model = 'BC')
- BeTSSi submarine (model = 'BCA')
- Mock shell (model = 'MS', model = 'Shirron', model = 'TAP')
- Barrel (model = 'Barrel')
- Torus (model = 'Torus')
- Cube (model = 'Cube')

## Overview of available methods
The main emphasis is acoustic scattering using plane wave (applyLoad = 'planeWave') in which the following cases are implemented
- Bistatic scattering (scatteringCase = 'BI')
- Monostatic scattering (scatteringCase = 'MS')
- Frequency sweep (scatteringCase = 'Sweep')


The following methods has been implemented (with available formulations)
- IGA using the infinite element method (method = 'IE')
	- The Bubnov--Galerkin Conjugated formulation (formulation = 'BGC')
	- The Petrov--Galerkin Conjugated formulation (formulation = 'PGC')
	- The Bubnov--Galerkin Unconjugated formulation (formulation = 'BGU')
	- The Petrov--Galerkin Unconjugated formulation (formulation = 'PGU')
- IGA using the IE method after Shirron (method = 'IENSG')
	- The Bubnov--Galerkin Conjugated formulation (formulation = 'BGC')
	- The Petrov--Galerkin Conjugated formulation (formulation = 'PGC')
	- The Bubnov--Galerkin Unconjugated formulation (formulation = 'BGU')
	- The Petrov--Galerkin Unconjugated formulation (formulation = 'PGU')
- IGA using absorbing boundary conditions (method = 'ABC')
	- Bayliss-Gunzburger-Turkel-operators (formulation = 'BGT')
- IGA using best approximation (method = 'BA'). Only available when analytic solution exist
	- Best approximation in L^2(\Gamma) (formulation = 'SL2E')
	- Best approximation in L^2(\Omega_a) (formulation = 'VL2E')
- IGA using the boundary element method (method = 'BEM')
	- Conventional boundary integral equation using collocation method (formulation = 'CCBIE')
	- Hyper-singular boundary integral equation using collocation method (formulation = 'CHBIE')
	- Burton-Miller formulation using collocation method (formulation = 'CBM')
	- Three regularized conventional boundary integral equations using collocation method (formulation = 'CCBIE1', formulation = 'CCBIE2' and formulation = 'CCBIE3')
  	
  	The Galerkin method may be used instead of the collocation method by replacing the prepended letter 'C' with 'G' (i.e. formulation = 'GCBIE'). Moreover, by appending the letter 'C' the CHIEF method is applied in an attempt to remove fictitious eigenfrequencies.
- IGA using Kirchhoff approximations (method = 'KDT')
	- A first order multiple bounce implementation (formulation = 'MS1')
	- (Not completed!) A second order multiple bounce implementation (formulation = 'MS2')
- The method of fundamental solutions (method = 'MFS')
	- Using fundamental solutions (Phi_k) as basis functions (formulation = 'PS')
	- Using spherical harmonics as basis functions (formulation = 'SS')
- IGA using Beam tracing (method = 'RT')


Instead of IGA (coreMethod = 'IGA') the following alternatives are implemented
- Bilinear isoparametric FEM (coreMethod = 'linear_FEM')
- Subparametric FEM with bilinear representation of geometry (coremethod = 'h_FEM')
- Isoparametric FEM (coreMethod = 'hp_FEM')

## Boundary conditions
The following boundary conditions are implemented
- Sound hard boundary condition (BC = 'SHBC')
- Sound soft boundary condition (BC = 'SSBC')
- Neumann-Neumann boundary condition (BC = 'NNBC')
- Neumann boundary condition (BC = 'NBC')
where the latter option is used for simulating manufactured solutions

## Additional parameters
In the getDefaultTaskValus.m file the addional parameters are described and set to some default values:

```Matlab
%% General settings
scatteringCase = 'BI';      % Bistatic scattering
applyLoad = 'planeWave';    % Acoustic scattering of plane wave
model = 'S1';               % Spherical shell of radius 1
method = 'IE';              % Method for handling the unbounded domain problem
coreMethod = 'IGA';         % Solution space
BC = 'SHBC';                % Boundary condition
loopParameters = {'M'};     % parameter study arr to be investigated
runTasksInParallel = false; % Run tasks in parallel
computeCondNumber = false;  % Compute the condition number of the global matrix
useROM = false;             % Use reduced order modeling
N_max = inf;                % number of terms in analytic solution for scattering on spherical shell

%% Storage settings
storeSolution               = false;    % Store the solution vector
storeFullVarCol             = false;    % Store all variable in the varCol variable collector
clearGlobalMatrices         = true;     % Clear memory consuming matrices

%% Mesh settings
M                   = 1;	% Mesh number
initMeshFactXi      = 1;	% initial number of knots in xi direction
initMeshFactEta 	= 1;	% initial number of knots in eta direction
initMeshFactZeta    = 1;	% initial number of knots in zeta direction
extraGP = 0;                % extra quadrature points

%% Settings for pre plotting (geometry and mesh visualisation)
plot2Dgeometry              = false;    % Plot cross section of mesh and geometry
plot3Dgeometry              = false;    % Plot visualization of mesh and geometry in 3D

%% Error computations
calculateSurfaceError       = false;	% Only if scatteringCase == 'Bi'
calculateSurfEnrgErr        = false;	% Only if scatteringCase == 'Bi'
calculateVolumeError        = false;	% Only if scatteringCase == 'Bi'

%% Settings for post plotting 
plotResultsInParaview       = false;	% Only if scatteringCase == 'Bi'
plotTimeOscillation         = false;	% Create 30 paraview files in order to visualize a dynamic result
plotMesh                    = false;	% Create additional Paraview files to visualize IGA mesh
plotFarField                = true;     % If false, plots the near-field instead
calculateFarFieldPattern    = true;     % Calculate far field pattern
alpha_s = NaN;                          % Aspect angle of incident wave
beta_s  = NaN;                          % Elevation angle of incident wave
alpha   = (0:0.5:360)*pi/180;           % Aspect angles of observation points
beta    = 0;                        	% Elevation angle of observation points
r       = 1;                            % radii for near-field evaluation. Assume by default plotFarField = true

%% Settings for the BEM (boundary element method)
useNeumanProj               = false;    % In BEM; project Neumann boundary conditions onto solution space
solveForPtot                = false;    % In BEM: solve for the total pressure (as opposed to the scattered pressure)
exteriorProblem             = true;     % In BEM: solve for the exterior problem (as opposed to the interior problem)
plotResidualError           = false;    % In Galerkin BEM; plot residual error of the BIE
extraGPBEM                  = 50;        	% extra quadrature points around singularities for BEM formulations
agpBEM                      = 1.4;       	% parameter for adaptiv Gauss point integration around singularities for BEM formulations
quadMethodBEM               = 'Adaptive2';  % In BEM: Quadrature method for handling weakly singular integrals
colBEM_C0                   = 0;        	% In collocation BEM: the scaling factor moving collocation points away from C0-lines
colMethod                   = 'Grev';     	% In collocation BEM: Location of collocation points
```

## Performing a study
A study may be perform over several sets of parameters. For example, using the following commands in the task script
```Matlab
M = [1,2,5];
degree = [2,5];
loopParameters = {'M','degree'};
```
will run separate simulations through all combinations of the mesh number 'M' and the polynomial degree 'degree'

## Disclaimer
No efforts has as of now been made to make the program user-friendly, and certain combinations of parameters/methods are not implemented or may be erroneous. Please report any bug to jonvegard89@yahoo.com
