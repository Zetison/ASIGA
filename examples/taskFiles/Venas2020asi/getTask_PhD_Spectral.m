function studies = getTask_PhD_Spectral()
% This study is based on Venas2020asi (Figure 7)
% Hetmaniuk2012raa is available at https://doi.org/10.1002/nme.4271

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% Set parameters
misc.scatteringCase = 'BI';

misc.model = 'S1';  % Spherical shell

misc.coreMethod = {'IGA'};

varCol = setS1Parameters('double',1);
msh.meshFile = 'createNURBSmesh_EL';
f = 1e3;
misc.omega = 2*pi*f;

ffp.alpha_s = 0;
ffp.beta_s = 0; 
ffp.calculateFarFieldPattern = 0;  

msh.M = 1:6;
msh.parm = 2;
msh.degree = 4;

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
% calculateVolumeError = 1;
err.calculateSurfaceError = 1;
misc.computeCondNumber = false;
misc.applyLoad = 'planeWave';

loopParameters = {'msh.M','msh.degree','misc.formulation','misc.coreMethod','misc.method'};

postPlot(1).xname        	= 'dofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).legendEntries 	= {};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
misc.formulation = {'SL2E'};
msh.M = 1:6; % 1:6
collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
misc.method = {'MFS'};
misc.formulation = {'PS'};
msh.M = 1:4; % 1:4
mfs.delta = 0.5;
% collectIntoTasks
