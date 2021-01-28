# ASIGA
The ASIGA (Acoustic Scattering using IsoGeometric Analysis) toolbox provides a framework for simulating acoustic scattering problems using isogeometric analysis (IGA). The toolbox contains a range of different method for solving acoustic scattering problems setting the ground work for research in this area. A range of geometries are included as well, in particular the BeTSSi submarine illustrated by the figure below. 
![Missing BeTSSi submarine nearfield plot](https://github.com/Zetison/ASIGA/blob/master/miscellaneous/BCA_p_tot.png?raw=true)
The main focus is to compute the so-called target strength (TS) for a given geometry as illustrated in the next figure.
![Missing target strength plot for the BeTSSi submarine](https://github.com/Zetison/ASIGA/blob/master/miscellaneous/BC_HWBC_MS_AS_E0_F1.png?raw=true)
The underlying governing equations and results produced by this toolbox can be found in \[[1](https://doi.org/10.1016/j.cma.2018.02.015), [2](https://www.sciencedirect.com/science/article/pii/S0045782519305559)\].

## Installation
Get the toolbox from GitHub (the e3Dss repository is needed for exact solutions and export_fig are needed for improved graphics for publications):
```
git clone https://github.com/altmany/export_fig
git clone https://github.com/Zetison/e3Dss 
git clone https://github.com/Zetison/ASIGA
```
In order to convert .iges-files to g2-files (which can be read using ASIGA), one needs the iges2go converter program from GoTools (which can be installed from https://github.com/SINTEF-Geometry/GoTools).
```
git clone --recursive https://github.com/SINTEF-Geometry/GoTools.git
git remote add Zetison https://github.com/Zetison/GoTools
git pull Zetison master

git clone git@github.com:Zetison/GoTools.git
cd GoTools
mkdir build
cd build
ccmake .. # Follow instructions
make
sudo make install
sudo ln -s $HOME/kode/GoTools/build/igeslib/app/iges2go /usr/local/bin/iges2go
iges2go < ~/kode/ASIGA/miscellaneous/FreeCADsphere.iges > ~/kode/ASIGA/miscellaneous/FreeCADsphere.g2

```

## Run program
The following will run the input script scriptName.m (given that studyName = 'scriptName'; is set in availableStudies.m)
```
matlab -nodisplay -nodesktop -nosplash -r "cd ASIGA, main; quit"
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
scatteringCase     = 'BI';        % Bistatic scattering
applyLoad          = 'planeWave'; % Set load. I.e.: 'planeWave', 'radialPulsation', 'pointPulsation', 'SimpsonTorus'
model              = 'S1';        % Spherical shell of radius 1
method             = 'IE';        % Method for handling the unbounded domain problem
coreMethod         = 'IGA';       % Solution space
BC                 = 'SHBC';      % Boundary condition
formulation        = 'BGU';       % Formulation
loopParameters     = {'M'};       % parameter study arr to be investigated
runTasksInParallel = false;       % Run tasks in parallel
computeCondNumber  = false;       % Compute the condition number of the global matrix
N_max              = inf;         % number of terms in analytic solution for scattering on spherical shell
progressBars       = false;        % Show progress bars for building system matrices
subFolderName      = '';          % sub folder in folder <folderName> in which results are stored

%% Storage settings
storeSolution       = false;    % Store the solution vector
storeFullVarCol     = false;    % Store all variable in the varCol variable collector
clearGlobalMatrices = true;     % Clear memory consuming matrices

%% Mesh settings
M                   = 1;	% Mesh number
initMeshFactXi      = 1;	% initial number of knots in xi direction
initMeshFactEta 	= 1;	% initial number of knots in eta direction
initMeshFactZeta    = 1;	% initial number of knots in zeta direction
extraGP             = 0;    % extra quadrature points
parm                = 1;    % Toggle different parameterizations of a geometric model

%% Settings for pre plotting (geometry and mesh visualisation)
prePlot.plot2Dgeometry      = false;       % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry      = false;       % Plot visualization of mesh and geometry in 3D
prePlot.storeFig            = false;       % Store pre plotted figure
prePlot.abortAfterPlotting  = false;       % Abort simulation after pre plotting
prePlot.view                = getView;     % Set view angle [azimuth,elevation]
prePlot.export_fig_name2D   = '';          % Name of exported figure using export_fig for 2D plots
prePlot.export_fig_name3D   = '';          % Name of exported figure using export_fig for 3D plots
prePlot.useCamlight         = true;        % Toggle camlight on

prePlot.plotControlPolygon  = false;       % Plot the control polygon for the NURBS mesh
prePlot.plotNormalVectors   = false;       % Plot the normal vectors for the NURBS mesh
prePlot.resolution       	  = [20,20,20];  % Number of evaluation points in the visualization for each element for each parametric direction
prePlot.color               = getColor(1); % Color for the NURBS surfaces
prePlot.plotAt              = true(3,2);   % For solids: toggle which surfaces to visualize
prePlot.alphaValue          = 1;           % Transparency value of NURBS surfaces
prePlot.colorFun            = NaN;         % Color the NURBS surfaces by the values of a colorFun function \R^d -> \R
prePlot.lineColor           = 'black';     % Mesh line color
prePlot.colorControlPolygon = 'red';       % Control polygon line color
prePlot.markerEdgeColor     = 'black';     % Control polygon edge marker color
prePlot.markerColor         = 'black';     % Control polygon marker color
prePlot.LineWidth           = 0.5;         % Width of lines
prePlot.elementBasedSamples = false;       % If true, sampling is based on distance rather than elements
prePlot.samplingDistance  	= NaN;         % Set sampling distance if elementBasedSamples = true
prePlot.title               = '';          % Set figure title
prePlot.axis                = 'off';       % Set axis() property
prePlot.xlabel              = 'x';         % Set x-axis label
prePlot.ylabel              = 'y';         % Set y-axis label
prePlot.zlabel              = 'z';         % Set z-axis label

%% Error computations
calculateSurfaceError       = false;	% Only if scatteringCase == 'Bi'
calculateSurfEnrgErr        = false;	% Only if scatteringCase == 'Bi'
calculateVolumeError        = false;	% Only if scatteringCase == 'Bi'
LpOrder                     = 2;        % Sets p for the L^p-norm

%% Settings for 1D far field evaluations
plotFarField                = true;     % If false, plots the near field instead
calculateFarFieldPattern    = true;     % Calculate far field pattern
alpha_s = NaN;                          % Aspect angle of incident wave
beta_s  = NaN;                          % Elevation angle of incident wave
alpha   = (0:0.5:360)*pi/180;           % Aspect angles of observation points
beta    = 0;                        	% Elevation angle of observation points
r       = 1;                            % radii for near-field evaluation. Assume by default plotFarField = true

%% Settings for post plotting 
postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= NaN;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).addCommands   	= [];

%% Settings for paraview
para.name                    = '';
para.plotResultsInParaview	 = false;	% Only if scatteringCase == 'Bi'
para.plotMesh              	 = true;	% Create additional Paraview files to visualize IGA mesh
para.plotP_inc               = true;
para.plotScalarField         = true;
para.plotTotField            = true; 
para.plotTotFieldAbs         = true; 
para.plotAnalytic            = true; 
para.plotTimeOscillation     = false;
para.plotDisplacementVectors = true;
para.computeGrad             = true;
para.plotError               = true; 
para.plotArtificialBoundary  = true;
para.extraXiPts              = 'round(20/2^(M-1))';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element

%% Settings for the BEM (boundary element method)
useNeumanProj       = false;    % In BEM; project Neumann boundary conditions onto solution space
solveForPtot     	= false;    % In BEM: solve for the total pressure (as opposed to the scattered pressure)
exteriorProblem  	= true;     % In BEM: solve for the exterior problem (as opposed to the interior problem)
plotResidualError	= false;    % In Galerkin BEM; plot residual error of the BIE
extraGPBEM         	= 50;        	% extra quadrature points around singularities for BEM formulations
agpBEM              = 1.4;       	% parameter for adaptiv Gauss point integration around singularities for BEM formulations
quadMethodBEM   	= 'Adaptive2';  % In BEM: Quadrature method for handling weakly singular integrals
colBEM_C0        	= 0;        	% In collocation BEM: the scaling factor moving collocation points away from C0-lines
colMethod         	= 'Grev';     	% In collocation BEM: Location of collocation points
internalPts      	= zeros(1,3);   % Internal points for the CHIEF method (Combined Helmholtz integral equation formulation)

%% Settings for the IEM (infinite element method)
r_a   	= NaN;    % Radial value for the artificial boundary
N     	= 3;        % Number of basis function in the radial direction for the IEM
IEbasis	= 'Chebyshev';

%% Settings for the MFS (method of fundamental solution)
delta = 0.1;    % Distance from the boundary to the internal source points

%% Settings for ROM (reduced order modelling)
useROM      = false;    % Toggle the usage of ROM
noVecsArr 	= 8;        % Number of derivatives at each point (including the 0th derivative
k_ROM    	= 3;        % Points at which to compute derivatives
useDGP      = true;     % Use a Derivative based Galerkin Projection (instead of interpolating techniques)
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
Certain combinations of parameters/methods are not implemented or may be erroneous. Please report any bug to jonvegard.venas@sintef.no
