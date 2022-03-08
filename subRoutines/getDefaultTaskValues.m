%% Default task values
% This script sets the default task values for each study

loopParameters     = {'msh.M'};   % parameter study array to be investigated
runTasksInParallel = false;       % Run tasks in parallel
subFolderName      = '';          % sub folder in folder <folderName> in which results are stored
saveStudies        = true;       % save ASIGA-struct into a .mat file
noCoresToUse       = Inf;         % Number of processors for parallel computations (Inf uses all available cores)
connectedParameters = {{}};       % Define set of loop parameters to be connected (i.e. connectedParameters = {{'msh.M','iem.N'}} assumes the arrays msh.M and iem.N to be of same size and loops through the elements in pairs)


%% Miscellaneous settings
misc.scatteringCase      = 'BI';        % Bistatic scattering
misc.applyLoad           = 'planeWave'; % Set load. I.e.: 'planeWave', 'radialPulsation', 'pointPulsation', 'SimpsonTorus'
misc.model               = 'S1';        % Spherical shell of radius 1
misc.method              = 'IE';        % Method for handling the unbounded domain problem
misc.coreMethod          = 'IGA';       % Solution space: 
misc.BC                  = 'SHBC';      % Boundary condition
misc.formulation         = 'BGU';       % Formulation
misc.computeCondNumber   = false;       % Compute the condition number of the global matrix
misc.N_max               = inf;         % number of terms in analytic solution for scattering on spherical shell
misc.progressBars        = false;       % Show progress bars for building system matrices
misc.solveForPtot     	 = false;       % In BEM or BA: solve for the total pressure (as opposed to the scattered pressure)
misc.plotResidualError	 = false;       % In Galerkin BEM; plot residual error of the BIE
misc.extraGP             = zeros(1,3);  % extra quadrature points
misc.r_a   	             = NaN;         % Radial value for the artificial boundary
misc.storeSolution       = false;       % Store the solution vector
misc.storeFullVarCol     = false;       % Store all variable in the varCol variable collector
misc.clearGlobalMatrices = true;        % Clear memory consuming matrices
misc.P_inc               = 1;           % Amplitude of incident wave
misc.exteriorProblem  	 = true;        % Solve for the exterior problem (as opposed to the interior problem)
misc.checkNURBSweightsCompatibility = true; % Check if the NURBS weights are compatible across patch interfaces

%% Mesh settings
msh.M                   = 1;	   % Mesh number
msh.initMeshFactXi      = 1;	   % initial number of knots in xi direction
msh.initMeshFactEta 	= 1;	   % initial number of knots in eta direction
msh.initMeshFactZeta    = 1;	   % initial number of knots in zeta direction
msh.parm                = 1;       % Toggle different parameterizations of a geometric model
msh.refineThetaOnly     = 0;       % For the ellipsoidal/spherical geometries, refine in the theta direction only
msh.degree              = 2;       % NURBS polynomial degree
msh.explodeNURBS        = false;   % Create patches from all C^0 interfaces
msh.x_0                 = [0,0,0]; % Translate center of model
msh.Xi                  = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
msh.nonLinearParam      = false;   % Use a non-linear parametrization to have equidistant control points (only implemented in the zeta-dir, and for geometrically linear parametrization)

%% Settings for pre plotting (geometry and mesh visualisation)
prePlot.plotFullDomain      = true;        % Plot volumetric domains
prePlot.plotSubsets         = {};          % Plot (surface) subsets (i.e. the artificial boundary Gamma_a) in paraview 
                                           % (examples include: 'Gamma','Gamma_a','yz','xz','xy','innerCoupling','outerCoupling','outer','inner','homDirichlet')
prePlot.plot2Dgeometry      = false;       % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry      = false;       % Plot visualization of mesh and geometry in 3D
prePlot.storeFig            = false;       % Store pre plotted figure
prePlot.abortAfterPlotting  = false;       % Abort simulation after pre plotting
prePlot.view                = getView;     % Set view angle [azimuth,elevation]
prePlot.export_fig_name2D   = '';          % Name of exported figure using export_fig for 2D plots
prePlot.export_fig_name3D   = '';          % Name of exported figure using export_fig for 3D plots
prePlot.useCamlight         = true;        % Toggle camlight on
preplot.camproj             = 'perspective'; % 'perspective' or 'orthographic'

plotParmDir = 1;
prePlot.plotControlPolygon  = false;       % Plot the control polygon for the NURBS mesh
prePlot.plotNormalVectors   = false;       % Plot the normal vectors for the NURBS mesh
prePlot.plotParmDir         = false;       % Plot arrows indication parametric directions
prePlot.plotGeometryInfo    = false;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
prePlot.plotArtificialBndry = true;        % Plot the artificial boundary for the IENSG method
prePlot.resolution       	= [20,20,20];  % Number of evaluation points in the visualization for each element for each parametric direction
prePlot.plotAt              = true(3,2);   % For solids: toggle which surfaces to visualize
prePlot.alphaValue          = 1;           % Transparency value of NURBS surfaces
prePlot.colorFun            = NaN;         % Color the NURBS surfaces by the values of a colorFun function \R^d -> \R
prePlot.lineColor           = 'black';     % Mesh line color
prePlot.colorControlPolygon = 'red';       % Control polygon line color
prePlot.markerEdgeColor     = 'black';     % Control polygon edge marker color
prePlot.markerColor         = 'black';     % Control polygon marker color
prePlot.LineWidth           = 0.1;         % Width of lines
prePlot.elementBasedSamples = false;       % If true, sampling is based on distance rather than elements
prePlot.samplingDistance  	= NaN;         % Set sampling distance if elementBasedSamples = true
prePlot.title               = '';          % Set figure title
prePlot.axis                = 'off';       % Set axis() property
prePlot.xlabel              = 'x';         % Set x-axis label
prePlot.ylabel              = 'y';         % Set y-axis label
prePlot.zlabel              = 'z';         % Set z-axis label
prePlot.pngResolution       = '-r200';     % Resolution of exported png images
prePlot.color               = [];          % Nx3 array with colors to use

%% Solver settings
sol.solver = 'LU';              % 'LU', 'gmres', 'cgs', 'bicgstab', 'bicgstabl', 'lsqr', 'bicg'
sol.preconditioner = 'diag';	% 'ilu', 'SSOR', 'diag'

%% Error computations
err.calculateSurfaceError = false;	% Only if scatteringCase == 'Bi'
err.calculateSurfEnrgErr  = false;	% Only if scatteringCase == 'Bi'
err.calculateVolumeError  = false;	% Only if scatteringCase == 'Bi'
err.LpOrder               = 2;      % Sets p for the L^p-norm

%% Settings for 1D far field evaluations
ffp.plotFarField                 = true;     % If false, plots the near field instead
ffp.calculateFarFieldPattern     = true;     % Calculate far field pattern
ffp.splineBasedNFPcalc           = false;    % Calculate near field pattern (NFP) directly from the splines
ffp.farFieldNormalPressFromSolid = true;
ffp.alpha_s = NaN;                          % Aspect angle of incident wave
ffp.beta_s  = NaN;                          % Elevation angle of incident wave
ffp.alpha   = (0:0.5:360)*pi/180;           % Aspect angles of observation points
ffp.beta    = 0;                        	% Elevation angle of observation points
ffp.r       = 1;                            % radii for near-field evaluation. Assume by default plotFarField = true
ffp.extraGP = [0,0,0];                      % Extra Gauss points used for the integration routine

%% Settings for post plotting 
postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';    % Examples include: 'p_Re', 'p_Im', 'abs_p', 'TS', 'error_pAbs', 'error_p', 'surfaceError', 'energyError', 'L2Error', 'H1Error', 'H1sError'
postPlot(1).plotResults  	= false;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= NaN;
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).addCommands   	= [];

%% Settings for paraview
para.name                    = '';
para.plotResultsInParaview	 = false;
para.plotMesh              	 = true;	% Create additional Paraview files to visualize IGA mesh
para.plotP_inc               = true;
para.plotScalarField         = true;
para.plotTotField            = true; 
para.plotTotFieldAbs         = true; 
para.plotAnalytic            = true; 
para.plotTimeOscillation     = false;
para.computeGrad             = true;
para.plotError               = true; 
para.plotSubsets             = {'Gamma','Gamma_a'}; % Plot (surface) subsets (i.e. the artificial boundary Gamma_a) in paraview 
                                                    % (examples include: 'Gamma','Gamma_a','yz','xz','xy','innerCoupling','outerCoupling','outer','inner','homDirichlet')
para.plotFullDomain          = true;
para.plotDisplacementVectors = true;
para.plotVonMisesStress      = true;
para.plotStressXX            = false;
para.plotStressYY            = false;
para.plotStressZZ            = false;
para.plotStressYZ            = false;
para.plotStressXZ            = false;
para.plotStressXY            = false;
para.i_MS                    = 1;   % Visualize of UU(:,i_MS) in paraview

para.extraXiPts              = 'round(20/2^(M-1))';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element

%% Settings for the BEM (boundary element method)
bem.useNeumanProj       = false;        % In BEM; project Neumann boundary conditions onto solution space
bem.extraGPBEM         	= 50;        	% extra quadrature points around singularities for BEM formulations
bem.agpBEM              = 1.4;       	% parameter for adaptiv Gauss point integration around singularities for BEM formulations
bem.quadMethodBEM   	= 'Adaptive2';  % In BEM: Quadrature method for handling weakly singular integrals
bem.colBEM_C0        	= 0;        	% In collocation BEM: the scaling factor moving collocation points away from C0-lines
bem.colMethod         	= 'Grev';     	% In collocation BEM: Location of collocation points
bem.internalPts      	= zeros(1,3);   % Internal points for the CHIEF method (Combined Helmholtz integral equation formulation)

%% Settings for the IEM (infinite element method)
iem.N     	 = 3;            % Number of basis function in the radial direction for the IEM
iem.IEbasis	 = 'Chebyshev';  % Choose between 'Chebyshev', 'Bernstein' and 'Lagrange'
iem.ie_Zeta  = [];           % Node/knot placements for basis function with local support in the infinite elements
iem.IElocSup = false;        % Toggle usage of radial shape functions in IE with local support
iem.p_ie     = NaN;          % Set polynomial order for radial shape functions
iem.s_ie     = NaN;          % Distrubution order for radial shape functions
iem.x_0      = zeros(1,3);   % The center of the prolate coordinate system of the infinite elemenets
iem.A_2      = eye(3);       % Rotation matrix for the prolate coordinate system 
iem.Upsilon  = 0;            % Parameter for prolate spheroidal coordinate system
iem.boundaryMethod = true;   % Attach infinite elements directly onto the scatterer for the IENSG formulation

%% Settings for the PML (perfectly matched layers)
pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
pml.sigmaType = 3;   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = gamma*xi^n, sigmaType = 3: sigma(xi) = gamma/(1-xi)^n, sigmaType = 4: sigma(xi) = gamma*(1/(1-xi)^n - 1), sigmaType = 5: sigma from Mi2021ilc: sigma(xi) = gamma*xi^n
pml.gamma = NaN;        % If ~isnan(pml.gamma) the matrices will be frequency independent (needed if useROM)
pml.t = NaN;         	% thickness of PML
pml.n = 1;            	% polynomial order
pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
pml.alpha = 30;         % constant used for the stretching function sigma(xi) for sigmaType = 5

%% Settings for the MFS (method of fundamental solution)
mfs.delta = 0.1;            % Distance from the boundary to the internal source points

%% Settings for RT (Ray tracing)
rt.N = 4; % Number of rays is approximately round(10^(rt.N/2))

%% Settings for ROM (reduced order modelling)
rom.useROM       = false;    % Toggle the usage of ROM
rom.noVecsArr 	 = 8;        % Number of derivatives at each point (including the 0th derivative
rom.k_ROM    	 = 3;        % Points at which to compute derivatives
rom.basisROMcell = {'DGP'};  % Basis for ROM ('Pade','Taylor','DGP','Hermite' or 'Bernstein')





