%% Default task values
% This script sets the default task values for each study

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
progressBars       = false;       % Show progress bars for building system matrices
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
refineThetaOnly     = 0;    % For the ellipsoidal/spherical geometries, refine in the theta direction only

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
prePlot.plotArtificialBndry = true;        % Plot the artificial boundary for the IENSG method
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


%% Solver settings
solver = 'LU';	        % 'LU', 'gmres', 'cgs', 'bicgstab', 'bicgstabl', 'lsqr', 'bicg'
preconditioner = 'ilu';	% 'ilu', 'SSOR'

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
r_a   	 = NaN;          % Radial value for the artificial boundary
N     	 = 3;            % Number of basis function in the radial direction for the IEM
IEbasis	 = 'Chebyshev';  % Choose between 'Chebyshev', 'Bernstein' and 'Lagrange'
ie_Zeta  = [];           % Node/knot placements for basis function with local support in the infinite elements
IElocSup = false;        % Toggle usage of radial shape functions in IE with local support
p_ie     = NaN;          % Set polynomial order for radial shape functions
s_ie     = NaN;          % Distrubution order for radial shape functions

%% Settings for the MFS (method of fundamental solution)
delta = 0.1;            % Distance from the boundary to the internal source points

%% Settings for ROM (reduced order modelling)
useROM      = false;    % Toggle the usage of ROM
noVecsArr 	= 8;        % Number of derivatives at each point (including the 0th derivative
k_ROM    	= 3;        % Points at which to compute derivatives
useDGP      = true;     % Use a Derivative based Galerkin Projection (instead of interpolating techniques)





