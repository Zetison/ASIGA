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

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha = (0:0.1:360)*pi/180;
beta = 30*pi/180;
f = 1000;

 
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 0;
prePlot.colorFun = @(v) abs(norm2(v)-1);
prePlot.resolution = [200,400];

postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= '';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

degree = 3; % FreeCAD exports to p = 3
M = 1; % 5
solveForPtot = true;
loopParameters = {'M','degree'};
collectIntoTasks

