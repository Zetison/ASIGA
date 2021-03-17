

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'S1';

misc.method = 'IE';
formulation = 'PGU';
misc.coreMethod = {'SEM','IGA'};
% misc.coreMethod = {'SEM'};
% misc.coreMethod = {'IGA'};
% f = [1e3,1e4];             % Frequency
f = 1e3;             % Frequency
% f = 3.98e4;             % Frequency

varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_EL';

M = 1; %3:6
N = [3,9];
N = 3;
N_max = N-1;
% alpha_s = 0;
% beta_s = pi/2;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
% alpha_s = pi;
% alpha_s = 0;
% beta_s = 0;
alpha = (0:0.5:360)*pi/180;
% alpha = 0;
% alpha = [0,90,180,270,360]*pi/180;
% beta = 30*pi/180;
beta = beta_s;
calculateVolumeError = 1;
err.err.calculateSurfaceError = 0;
computeCondNumber = true;
clearGlobalMatrices = false;
% misc.applyLoad = 'pointPulsation'; % with analytic solution for arbitrary geometries

postPlot(1).xname        	= 'dofsAlg';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'degree';
postPlot(1).legendEntries 	= {'misc.coreMethod','f'};
postPlot(1).subFolderName 	= NaN;
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(2) = postPlot(1);
postPlot(2).axisType = 'loglog';
postPlot(2).xname = 'dofs';
postPlot(2).yname = 'cond_number';

postPlot(3) = postPlot(2);
postPlot(3).yname = 'energyError';
postPlot(3).xname = 'tot_time';

postPlot(4) = postPlot(3);
postPlot(4).xname = 'timeSolveSystem';

postPlot(5) = postPlot(4);
postPlot(5).xname = 'timeBuildSystem';
    

calculateFarFieldPattern = 0;
prePlot.plot3Dgeometry = 0;
degree = 4:19;
degree = 4:10;
parm = 2;

loopParameters = {'degree','N','misc.coreMethod','f'};
collectIntoTasks


