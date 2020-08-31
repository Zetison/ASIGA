function studies = getTask_M31()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M31';

method = {'BEM'};
formulation = {'GCBIE'};

varCol = setM31Parameters(3);
varCol{1}.meshFile = 'createNURBSmesh_M31';
f = 1e2;             % Frequency

M = 3:4;
% M = 1;
degree = 4;
alpha = (0:0.5:360)*pi/180;
beta = 0;
alpha_s = 240*pi/180;
beta_s = 0;
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};

prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'on';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting  = 0;                % Abort simulation after pre plotting

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

collectIntoTasks


