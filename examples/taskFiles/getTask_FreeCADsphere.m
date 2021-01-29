function studies = getTask_FreeCADsphere()
% This study imports .g2 files from FreeCAD (.g2 file obtained from a iges
% files, see README.md). Remember that ASIGA assumes the model to be in
% meters, so make sure to export as such (alternatively use the scaleNURBS()-function).

counter = 1;
studies = cell(0,1);
getDefaultTaskValues


%% IE simulation
scatteringCase = 'BI';
% BC = 'NNBC';
model = 'FreeCADsphere';  % Spherical shell

coreMethod = {'IGA'};
method = {'BEM'};
formulation = {'CCBIE'};
noDomains = 1;
varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_FreeCADsphere';

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha = (0:0.1:360)*pi/180;
beta = 30*pi/180;
f = 1000;

 
calculateSurfaceError = 1;
calculateFarFieldPattern = 1;
prePlot.abortAfterPlotting  = 1;       % Abort simulation after pre plotting
prePlot.plot3Dgeometry = 1;
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

postPlot(2) = postPlot(1);
postPlot(2).xname        	= 'alpha';
postPlot(2).yname        	= 'error_p';
postPlot(2).axisType      	= 'semilogy';

postPlot(3) = postPlot(1);
postPlot(3).noXLoopPrms = 1;
postPlot(3).xname        	= 'nepw';
postPlot(3).yname        	= 'surfaceError';
postPlot(3).xLoopName     	= 'M';
postPlot(3).axisType      	= 'loglog';

degree = 3; % FreeCAD exports to p = 3
M = 1:2; % 5
solveForPtot = true;
loopParameters = {'M','method'};
% collectIntoTasks

method = {'BA'};
formulation = {'SL2E'};
M = 1:4; % 5
collectIntoTasks

