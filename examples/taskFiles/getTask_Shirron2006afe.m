function studies = getTask_Shirron2006afe()
% This study is based on Shirron2006afe (Fig. 2)
% Shirron2006afe is available athttps://doi.org/10.1016/j.cma.2006.07.009

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'Shirron2006afe';  % Spherical shell
misc.coreMethod = {'C0_IGA'};
% misc.coreMethod = {'IGA'};
misc.applyLoad = 'planeWave';
% misc.applyLoad = 'radialPulsation';

% misc.method = {'IENSG'};
% misc.method = {'IE'};
misc.method = {'PML'};
misc.BC = {'SHBC'};

ffp.calculateFarFieldPattern    = true;     % Calculate far field pattern
ffp.alpha_s = 0;                            % Aspect angle of incident wave
ffp.beta_s  = -pi/2;                        % Elevation angle of incident wave
ffp.alpha   = 0;                            % Aspect angles of observation points
ffp.beta = linspace(-pi/2,pi/2,1000);   
ffp.r = 1;                            % radii for near-field evaluation.
ffp.paramPts = {zeros(1000,3), []};
ffp.paramPts{1}(:,2) = linspace(0,1,1000);
ffp.splineBasedNFPcalc = true;

a = 1; % Midsurface radius
varCol{1} = struct('media', 'fluid', ...
                  'R_i', a, ...
                  't',   a, ...
                  'c_f', 1524, ...
                  'rho', 1000);
iem.N = 4; % 9
msh.meshFile = 'createNURBSmesh_EL';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.refineThetaOnly = true;
if strcmp(misc.method{1},'PML')
    misc.formulation = {'GSB'};
%     misc.formulation = {'STD'};
    varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
    varCol{1}.refinement = @(M) [0, round(pi/0.5*1.5*10/2), M, M];
%     varCol{1}.refinement = @(M) [0, 0, 4, M];
else
    misc.formulation = {'BGU'};
    varCol{1}.refinement = @(M) [0, round(pi/0.5*1.5*10/2), M, M];
%     varCol{1}.refinement = @(M) [0, 0, 4];
end

msh.degree = 4;
msh.M = 4; % 4
misc.extraGP = [7,0,0];    % extra quadrature points
% extraGP = [7,7,0];    % extra quadrature points
warning('off','NURBS:weights')

k = 10/a;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.25*a;

pml.eps = 1e9*eps;      % choosing eps = eps yields machine precicion at Gamma_b, but requires more "radial" elements in the PML to resolve the rapid decay function
pml.sigmaType = [1,2];  % sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
pml.n = 2;
pml.gamma = 5;
pml.t = 0.25*a;         % thickness of PML
pml.dirichlet = 0;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

msh.parm = 1;
err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;

ffp.calculateFarFieldPattern = 1;


prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
misc.computeCondNumber = 0;

postPlot(1).xname           = 'M';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'msh.M';

postPlot(2) = postPlot(1);
postPlot(2).plotResults  	= true;
postPlot(2).printResults 	= 1;
postPlot(2).xScale       	= 180/pi;
postPlot(2).xname           = 'beta';
postPlot(2).yname           = 'abs_p';
% postPlot(2).yname           = 'TS';
postPlot(2).axisType      	= 'plot';
postPlot(2).noXLoopPrms   	= 0;

% postPlot(3) = postPlot(2);
% postPlot(3).axisType = 'semilogy';
% postPlot(3).yname = 'error_p';

loopParameters = {'msh.M','misc.method','pml.sigmaType'};
para.plotResultsInParaview	= 0;
para.extraXiPts              = '30';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(1/2^(M-1))';   % Extra visualization points in the zeta-direction per element
% collectIntoTasks

%% PML parameter study
msh.degree = 2;
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, 2^(M-4)-1, 2^(M-4)-1];
msh.M = 6:7; 
% msh.M = 2:3; 

k = [5,10]/a;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.25*a;

misc.progressBars = 0;
runTasksInParallel = 1;       % Run tasks in parallel
misc.checkNURBSweightsCompatibility = false;

pml.sigmaType = 2:4;  % sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n
pml.n = [1,2];
% pml.gamma = NaN;
pml.gamma = linspace2(0,5,100); % pml.t*k = 2.5
pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

ffp.calculateFarFieldPattern = 0;

postPlot(1).xname           = 'pml.gamma';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType      	= 'semilogy';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'pml.gamma';
postPlot(2) = [];
loopParameters = {'msh.M','misc.method','pml.gamma','pml.sigmaType','pml.n','misc.omega'};
collectIntoTasks

misc.method = {'IE'};
misc.formulation = {'BGU'};
pml.sigmaType = NaN;
pml.n = NaN;
pml.gamma = [pml.gamma(1),pml.gamma(end)];
loopParameters = {'msh.M','misc.method','pml.gamma'};
collectIntoTasks

misc.method = {'BA'};
misc.formulation = {'SL2E'};
collectIntoTasks
