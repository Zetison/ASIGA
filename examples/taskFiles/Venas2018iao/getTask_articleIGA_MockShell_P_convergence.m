function studies = getTask_articleIGA_MockShell_P_convergence()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 21-23 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)


counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.model = 'MS'; % Mock shell misc.model after Ihlenburg

% misc.method = {'IE','IENSG'};
% formulation = {'BGU','PGC','BGC'};
misc.method = 'IE';
formulation = {'BGU','PGU','PGC','BGC'};
formulation = {'BGU'};
% misc.method = 'BA';
% formulation = {'VL2E'};

varCol = setMockShellParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M3';
varCol{1}.L = pi/2;

IEbasis = {'Chebyshev','Bernstein','Lagrange'};
IEbasis = {'Bernstein'};
M = 1:6;
M = 1:5;
N = [1,2,3,6,9];
N = [1,2];
% N = 1;
degree = 2;
parm = 1;
runTasksInParallel = 0;
plotResultsInParaview = 0;
err.calculateSurfaceError = 1;	% Only for spherical shell and if misc.scatteringCase == 'Bi'
LpOrder = 2; % For error calculation in calcSurfError()

calculateVolumeError  = 1;	% Only for spherical shell and if misc.scatteringCase == 'Bi'
calculateFarFieldPattern = 0;
computeCondNumber = 1;
misc.applyLoad = 'pointPulsation';
BC = 'NBC';
warning('off','NURBS:weights')

k = 10;
c_f = 1500;

f = k*c_f/(2*pi);             % Frequency

prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry = 0;  % Plot cross section of mesh and geometry
prePlot.abortAfterPlotting = 0;  % Plot cross section of mesh and geometry
prePlot.plotControlPolygon = 0;  % Plot cross section of mesh and geometry


postPlot(1).xname       	= 'dofs';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).addCommands   	= [];


loopParameters = {'formulation','M', 'N', 'IEbasis'};
% loopParameters = {'M'};
collectIntoTasks
