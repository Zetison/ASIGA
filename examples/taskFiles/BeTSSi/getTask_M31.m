function studies = getTask_M31()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M31';

misc.method = {'BEM'};
misc.formulation = {'GCBIE'};

varCol = setM31Parameters(3);
varCol{1}.meshFile = 'createNURBSmesh_M31';
f = 1e2;             % Frequency
misc.omega = 2*pi*f;

msh.M = 3:4;
% msh.M = 1;
msh.degree = 4;
ffp.alpha = (0:0.5:360)*pi/180;
ffp.beta = 0;
ffp.alpha_s = 240*pi/180;
ffp.beta_s = 0;
misc.solveForPtot = true;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};

prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
% prePlot.axis = 'on';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting  = 1;                % Abort simulation after pre plotting

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','msh.parm','misc.formulation','.msh.M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

collectIntoTasks


