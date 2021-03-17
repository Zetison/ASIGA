function studies = getTask_Avery2007ffs()
% This study is based on Simpson2014aib and correspond to Figure 6 in Venas2018iao
% Simpson2014aib is available at https://doi.org/10.1016/j.cma.2013.10.026
% Venas2018iao   is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'SS'; % Simpson sphere
misc.method = {'IE'};
BC = 'SHBC';
formulation = {'PGU','PGC','BGU','BGC'};
% formulation = {'BGU'};
misc.coreMethod = 'IGA';
computeCondNumber = 0;
runTasksInParallel = 0;
progressBars = false;        % Show progress bars for building system matrices

varCol = setSSParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_EL';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 2;                 % Wave number for Simpson2014aib
omega = c_f*k;         % Angular frequency
f = omega/(2*pi);      % Frequency

M = 3;

alpha_s = pi;
beta_s = 0;  
alpha = (0:0.5:360)*pi/180;
prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [20,20,0];

plotFarField = false; 
r = 5; % radii for near-field evaluation
degree = 3;

N = 4;

loopParameters = {'M','misc.method','formulation'};
err.calculateSurfaceError = 1;
parm = 1;

postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= [];

postPlot(2) = postPlot(1);
postPlot(2).noXLoopPrms     = 0;
postPlot(2).legendEntries   = {'misc.method','M','formulation'};
postPlot(2).lineStyle       = '-';
postPlot(2).xname           = 'alpha';
postPlot(2).yname           = 'error_pAbs';
postPlot(2).axisType        = 'semilogy';
postPlot(2).xScale          = 180/pi;
postPlot(2).yScale          = 1/100;
postPlot(2).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

postPlot(3) = postPlot(2);
postPlot(3).yname       = 'abs_p';
postPlot(3).axisType	= 'plot';
postPlot(3).yScale      = 1;
postPlot(3).addCommands = [];

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
useNeumanProj = 0;
solveForPtot = 1;
% formulation = {'SL2E','VL2E'};
formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
M = 1:3;
% formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};
colBEM_C0 = 0;
solveForPtot = true;
collectIntoTasks