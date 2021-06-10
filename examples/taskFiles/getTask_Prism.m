function studies = getTask_Prism()
% This study is based on Shirron2006afe (Fig. 2)
% Shirron2006afe is available athttps://doi.org/10.1016/j.cma.2006.07.009

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
misc.scatteringCase = 'MS';
misc.model = 'Prism';  % Spherical shell
misc.coreMethod = {'C0_IGA'};
% misc.coreMethod = {'IGA'};
misc.formulation = {'GSB'};
% misc.applyLoad = 'planeWave';
misc.applyLoad = 'pointPulsation';

% misc.method = {'IENSG'};
misc.method = {'PML'};
misc.BC = {'SHBC'};
misc.BC = {'NBC'};

ffp.alpha = 0;                            % Aspect angle of incident wave
ffp.beta  = -pi/2;                        % Elevation angle of incident wave

a = 1;
varCol{1} = struct('media', 'fluid', ...
                  'L', [4,1,2], ...
                  't',   0.1, ...
                  't_fluid', 0.25*a, ...
                  'c_f', 1500, ...
                  'rho', 1000);
              
msh.meshFile = 'createNURBSmesh_Prism';
varCol{1}.refinement = @(M) [0,0,0,2^(M-1)-1];

msh.degree = 2;
msh.M = 1:5;
msh.M = 6;
warning('off','NURBS:weights')

k = 1/a;
lambda = 2*pi/k;
f = k*varCol{1}.c_f/(2*pi);
misc.omega = 2*pi*f;
misc.r_a = 1.25*a;

pml.sigmaType = 3;
pml.n = 1;
pml.t = 0.25*a;         % thickness of PML
% pml.t = 0.125*a;         % thickness of PML
pml.gamma = 1/(pml.t*k(1));
pml.dirichlet = 1;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)

err.calculateSurfaceError = 1;
err.calculateVolumeError  = 0;

ffp.calculateFarFieldPattern = 0;


prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 1;
prePlot.plot2Dgeometry = 0;
prePlot.plotGeometryInfo = 0;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
prePlot.plotParmDir      = 0;       % Plot arrows indication parametric directions
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
prePlot.resolution = [0,0,0];
misc.computeCondNumber = 0;


para.plotResultsInParaview	 = 1;
para.extraXiPts              = '0';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = '0';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = '0';   % Extra visualization points in the zeta-direction per element

postPlot(1).xname           = 'surfDofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '-*';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'msh.M';
loopParameters = {'msh.M','msh.degree','misc.method','misc.formulation'};
collectIntoTasks
% 
misc.method = {'BA'};
misc.formulation = {'VL2E'};
% collectIntoTasks

misc.method = {'BA'};
misc.formulation = {'SL2E'};
% collectIntoTasks


