function studies = getTask_Barrel_conv()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

saveStudies = false;

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'Barrel';
BC = 'NBC';
misc.method = {'BEM'};
formulation = {'CCBIE'};
misc.method = {'BA'};
formulation = {'SL2E'};
% formulation = {'CCBIE','GBM'};

misc.applyLoad          = 'pointPulsation'; % Acoustic scattering of plane wave
varCol = setBarrelParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_Barrel';
f = 1.5e2;             % Frequency

M = 1:6;
% M = 4;
degree = 2;
parm = 1:2;
alpha = (0:0.5:360)*pi/180;
beta = 0;

solveForPtot = false;
err.calculateSurfaceError       = 1;	% Only if misc.scatteringCase == 'Bi'

warning('off','NURBS:weights')
loopParameters = {'M','parm','f','misc.method','formulation'};

prePlot.plot3Dgeometry = 0;
% prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 0;                % Abort simulation after pre plotting

postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','formulation','M','N','IEbasis'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= [];

collectIntoTasks


