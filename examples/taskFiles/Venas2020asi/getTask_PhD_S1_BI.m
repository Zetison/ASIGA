misc.scatteringCase = 'BI';

misc.model = 'S1';  % Spherical shell

misc.coreMethod = 'IGA';
varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_EL';


f = 20e3;

alpha_s = 0;
beta_s = 0;   

M = 3;
postPlot(1).xname        	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'plot';
postPlot(1).lineStyle    	= '-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'N';
postPlot(1).legendEntries 	= {'misc.coreMethod','f'};
postPlot(1).subFolderName 	= NaN;
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;

postPlot(2) = postPlot(1);
postPlot(2).noXLoopPrms = 0;
postPlot(2).lineStyle = '-';
postPlot(2).xname = 'alpha';
postPlot(2).yname = 'error_pAbs';
postPlot(2).axisType = 'semilogy';
postPlot(2).xScale = 180/pi;

postPlot(3) = postPlot(2);
postPlot(3).yname = 'error_p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
formulation = '';
M = 3;
prePlot.plot3Dgeometry = 0;
degree = 2;
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
misc.applyLoad = 'planeWave';
parm = 1;
N = 3:4;

loopParameters = {'N','M','misc.method'};
collectIntoTasks
