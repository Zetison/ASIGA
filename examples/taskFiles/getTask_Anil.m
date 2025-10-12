function studies = getTask_Anil()
% This study is based on Shirron2006afe (Fig. 2)
% Shirron2006afe is available athttps://doi.org/10.1016/j.cma.2006.07.009

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

saveStudies = false;


%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'Anil';  % Spherical shell
% misc.coreMethod = {'C0_IGA'};
misc.coreMethod = {'IGA'};
% misc.formulation = {'GSB'};
misc.applyLoad = 'planeWave';
% misc.applyLoad = 'radialPulsation';

% misc.method = {'IENSG'};
% misc.method = {'IE'};
% misc.method = {'PML'};
ffp.plotFarField = 1;     % If false, plots the near field instead
% ffp.splineBasedNFPcalc = [0,1];
ffp.splineBasedNFPcalc = 0;
if 0
    misc.method = {'PML','NIL'};
    misc.formulation = {'GSB', 'SPH'};
    pml.sigmaType = [3,1];
    pml.dirichlet = [1,0];	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
else
    misc.method = {'NIL'};
    misc.formulation = {'SPH'};
    pml.sigmaType = 1;
    pml.dirichlet = 0;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
end
ffp.paramPts = {[], [linspace(0,1,361).', 0.5*ones(361,1), ones(361,1)]};
misc.BC = {'SHBC'};
connectedParameters = {{'misc.method','misc.formulation','pml.sigmaType','pml.dirichlet'}};

a = 1;
varCol{1} = struct('media', 'fluid', ...
                  'R', a, ...
                  't',   a, ...
                  'c_f', 1500, ...
                  'rho', 1000);
iem.N = 4; % 9
msh.meshFile = 'createNURBSmesh_EL';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.refineThetaOnly = 0;
% varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-4)-1, max(2^(M-4)-1,3)];
msh.autoRefine          = true;
% varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 31];

msh.degree = [2,3];
msh.degree = 4;
% msh.M = 5:7;
msh.M = 4;
msh.explodeNURBS = 0;   % Create patches from all C^0 interfaces


ffp.calculateFarFieldPattern     = true;     % Calculate far field pattern

ffp.beta = 0;
% ffp.alpha = (0:0.05:360)*pi/180;
ffp.alpha = (0:1:360)*pi/180;
ffp.alpha_s = pi;
ffp.beta_s = 0;

warning('off','NURBS:weights')

k = 4;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.4*a;

pml.n = 1;
pml.t = 0.6*a;         % thickness of PML
pml.gamma = 1/(pml.t*k(1));

msh.parm = 1;
err.calculateSurfaceError = 0;
err.calculateVolumeError  = 1;

prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
misc.computeCondNumber = 0;
prePlot.plotFullDomain   = 1;        % Plot volumetric domains
prePlot.view             = [0,90];
prePlot.plotSubsets      = {'xy'};

misc.progressBars = 0;
runTasksInParallel = 0;       % Run tasks in parallel
misc.checkNURBSweightsCompatibility = false;

pml.eps = eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
% pml.sigmaType = 1:4;  % sigmaType = 1: sigma(xi) = xi*(exp(gamma*xi)-1), sigmaType = 2: sigma(xi) = C*xi^n
% pml.n = [1,2];
% pml.t = linspace2(0,0.25*a,100);
% pml.t = 0.25*a;
% pml.gamma = NaN;
% pml.gamma = linspace2(0,10,100); % pml.t*k = 2.5

para.plotResultsInParaview	 = 0;

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'abs_p';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
% postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
% postPlot(1).subFolderName   = '~/Dropbox/Apps/Overleaf/createFigures/data/PMLarticle';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
% postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);
postPlot(2) = postPlot(1);

postPlot(2).axisType        = 'semilogy';
postPlot(2).yname        	= 'error_p';
postPlot(2).yScale          = 1/100;
loopParameters = {'misc.method'};
% collectIntoTasks

%% Run convergence analysis

msh.M = 1:7;
msh.refineThetaOnly = 1;
ffp.calculateFarFieldPattern = 0;
para.plotResultsInParaview	 = 0;
ffp.alpha_s = pi;
ffp.beta_s = pi/2;
msh.degree = 4;

prePlot.plot3Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
misc.computeCondNumber = 0;
prePlot.plotFullDomain   = 0;        % Plot volumetric domains
prePlot.view             = [0,0];
prePlot.plotSubsets      = {'xz'};

if 1
    misc.method = {'PML','NIL'};
    misc.formulation = {'GSB', 'SPH'};
    pml.sigmaType = [3,1];
    pml.dirichlet = [1,0];	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    loopParameters = {'misc.method','msh.M','msh.degree','misc.coreMethod','pml.sigmaType'};
else
%     misc.method = {'NIL'};
%     misc.formulation = {'SPH'};
%     pml.sigmaType = 1;
%     pml.dirichlet = 0;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    misc.method = {'PML'};
    misc.formulation = {'GSB'};
    pml.sigmaType = [1,2,3,4,5];
    pml.n = 1:3;
    pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
    connectedParameters = {{'misc.method','misc.formulation','pml.dirichlet'}};
    loopParameters = {'misc.method','msh.M','msh.degree','misc.coreMethod','pml.sigmaType','pml.n'};
end

if msh.refineThetaOnly
    ffp.extraGP = [59,0,0];
    misc.extraGP = [10,0,0];    % extra quadrature points
    % extraGP = [7,7,0];    % extra quadrature points
end

postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addSlopes       = false;
postPlot(2) = postPlot(1);
postPlot(2).xname       	= 'varCol{1}.nepw';
postPlot(2).addSlopes       = true;
collectIntoTasks

msh.M = 1:7;
misc.coreMethod = {'C0_IGA'};
collectIntoTasks
