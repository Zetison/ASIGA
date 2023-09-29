function studies = getTask_Venas2020asi_Figure21(M_0)
% This study is based on Venas2020asi (Figure 21)
% Venas2020asi is available at http://hdl.handle.net/11250/2640443
% Note that the PhD thesis has a typo in saying this simulation uses p=4
% when the true value was p=2

if nargin < 1
    M_0 = 5; % 5
end
counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% Set parameters
misc.scatteringCase = 'BI';

misc.model = 'S1';  % Spherical shell

misc.coreMethod = {'IGA'};

varCol = setS1Parameters('double',1);
msh.meshFile = 'createNURBSmesh_EL';
msh.autoRefine = true;

f = 1e3;
misc.omega = 2*pi*f;

ffp.alpha_s = 0;
ffp.beta_s = 0; 
msh.parm = 1;
msh.degree = 2;
ffp.calculateFarFieldPattern = 1;  

warning('off','NURBS:weights')

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
% calculateVolumeError = 1;
err.calculateSurfaceError = 1;
misc.computeCondNumber = false;
prePlot.abortAfterPlotting  = 1;       % Abort simulation after pre plotting
misc.applyLoad = 'planeWave';

loopParameters = {'msh.M','msh.degree','misc.formulation','misc.coreMethod','misc.method'};

postPlot(1).xname        	= 'surfDofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
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
msh.M = 1:M_0+1; % 1:6
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
misc.method = {'MFS'};
misc.formulation = {'PS'};
msh.M = 1:M_0; % 1:5
% mfs.delta = 0.32778;
mfs.delta = 0.5;
collectIntoTasks


