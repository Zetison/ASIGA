scatteringCase = 'BI';

model = 'S1';  % Spherical shell

coreMethod = {'IGA'};

varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_EL';
f = 1e3;

alpha_s = 0;
beta_s = 0;   

M = 1:6;
parm = 2;

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
degree = 4;
% calculateVolumeError = 1;
calculateSurfaceError = 1;
computeCondNumber = false;
calculateFarFieldPattern = 0;
applyLoad = 'planeWave';

loopParameters = {'M','degree','formulation','coreMethod','method'};

postPlot(1).xname        	= 'dofs';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {};
postPlot(1).subFolderName 	= NaN;
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
method = {'BA'};
formulation = {'SL2E'};
M = 1:6;
collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
formulation = {'PS'};
M = 1:4;
delta = 0.5;
collectIntoTasks


