function studies = getTask_Barrel()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'Barrel';
BC = 'SHBC';
method = {'BEM'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};

varCol = setBarrelParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_Barrel';
f = 1.5e3;             % Frequency

M = 5:6;
% M = 5;
degree = 2;
beta = 0;
parm = 1;
alpha = (0:0.5:360)*pi/180;

quadMethodBEM = 'Simpson';
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};

prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 0;                % Abort simulation after pre plotting

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


method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks


