function studies = getTask_FreeCAD()
% This study imports .g2 files from FreeCAD (.g2 file obtained from a iges
% files, see README.md). Remember that ASIGA assumes the model to be in
% meters, so make sure to export as such (alternatively use the scaleNURBS()-function).

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'BI';
% BC = 'NNBC';
model = 'FreeCAD';  % Spherical shell

coreMethod = {'IGA'};
method = {'BEM'};
formulation = {'CCBIE'};
noDomains = 1;
varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_FreeCAD';


applyLoad = 'planeWave';
BC = 'SHBC';
alpha_s = 240*pi/180;
beta_s = 0;
alpha = (0:0.1:360)*pi/180;
beta = 0;
f = 100;


calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = false;       % Abort simulation after pre plotting
prePlot.plotNormalVectors = true;
prePlot.plot3Dgeometry = 1;
prePlot.plotControlPolygon = 1;
prePlot.resolution = [100,100];

postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'polar';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale          = 180/pi;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

para.name                    = '';
para.plotResultsInParaview	 = 1;	% Only if scatteringCase == 'Bi'

degree = 3; % FreeCAD exports to p = 3
M = 2:3; % 5
solveForPtot = true;
loopParameters = {'M','method'};
% collectIntoTasks

method = {'KDT'};
formulation = {'MS1'};
solveForPtot = false;
para.plotResultsInParaview = 0;	% Only if scatteringCase == 'Bi'
collectIntoTasks