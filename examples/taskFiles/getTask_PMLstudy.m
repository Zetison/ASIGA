function studies = getTask_PMLstudy()
% This study is based on Shirron2006afe (Fig. 2)
% Shirron2006afe is available athttps://doi.org/10.1016/j.cma.2006.07.009

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
misc.scatteringCase = 'MS';
misc.model = 'PMLstudy';  % Spherical shell
misc.coreMethod = {'C0_IGA'};
% misc.coreMethod = {'IGA'};
misc.formulation = {'GSB'};
misc.applyLoad = 'planeWave';
% misc.applyLoad = 'radialPulsation';

% misc.method = {'IENSG'};
% misc.method = {'IE'};
misc.method = {'PML'};
misc.BC = {'SHBC'};

ffp.alpha = 0;                            % Aspect angle of incident wave
ffp.beta  = -pi/2;                        % Elevation angle of incident wave

a = 1;
varCol{1} = struct('media', 'fluid', ...
                  'R_i', a, ...
                  't',   a, ...
                  'c_f', 1524, ...
                  'rho', 1000);
iem.N = 4; % 9
msh.meshFile = 'createNURBSmesh_EL';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.refineThetaOnly = 1;
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, max(2^(M-4)-1,3)];
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 31];

msh.degree = [2,3];
msh.degree = 2;
msh.M = 5:7;
% msh.M = 5;
misc.extraGP = [7,0,0];    % extra quadrature points
% extraGP = [7,7,0];    % extra quadrature points
warning('off','NURBS:weights')

k = 10/a;
k = 1/a;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.25*a;

pml.sigmaType = 4;
pml.n = 1;
pml.gamma = 1/(pml.t*k);
pml.t = 0.25*a;         % thickness of PML
pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

msh.parm = 1;
err.calculateSurfaceError = 0;
err.calculateVolumeError  = 1;

ffp.calculateFarFieldPattern = 0;


prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
misc.computeCondNumber = 0;

postPlot(1).xname           = 'surfDofs';
postPlot(1).yname        	= 'L2Error';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-*';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'msh.M';
loopParameters = {'msh.M','msh.degree','misc.method'};
% collectIntoTasks

iem.N = 4;
misc.method = {'IE'};
misc.formulation = {'BGU'};
% collectIntoTasks

misc.method = {'BA'};
misc.formulation = {'VL2E'};
% collectIntoTasks

%% sigmaType parameter study
misc.method = {'PML'};
misc.formulation = {'GSB'};
msh.degree = 2;
msh.refineThetaOnly = true;
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
msh.M = 5:7; 
% msh.M = 1:2; 

k = [5,10]/a;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.25*a;

misc.progressBars = 0;
runTasksInParallel = 1;       % Run tasks in parallel
misc.checkNURBSweightsCompatibility = false;

pml.eps = eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
pml.sigmaType = 1:4;  % sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
pml.n = [1,2];
% pml.gamma = NaN;
pml.gamma = linspace2(0,10,100); % pml.t*k = 2.5
pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

ffp.calculateFarFieldPattern = 0;

postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
postPlot(2).xname           = 'pml.gamma';
postPlot(2).yname        	= 'L2Error';
postPlot(2).plotResults  	= 1;
postPlot(2).printResults 	= 1;
postPlot(2).axisType      	= 'semilogy';
postPlot(2).lineStyle    	= '-';
postPlot(2).xScale       	= 1;
postPlot(2).yScale       	= 1;
postPlot(2).fileDataHeaderX	= [];
postPlot(2).noXLoopPrms   	= 1;
postPlot(2).xLoopName     	= 'pml.gamma';
loopParameters = {'msh.M','misc.method','pml.gamma','pml.sigmaType','pml.n','misc.omega'};
collectIntoTasks

iem.N = 4;
misc.method = {'IE'};
misc.formulation = {'BGU'};
pml.sigmaType = NaN;
pml.n = NaN;
pml.gamma = [pml.gamma(1),pml.gamma(end)];
loopParameters = {'msh.M','misc.method','pml.gamma','misc.omega'};
collectIntoTasks

misc.method = {'BA'};
misc.formulation = {'VL2E'};
collectIntoTasks
