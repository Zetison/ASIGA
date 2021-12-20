function studies = getTask_Barrel()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.applyLoad = 'pointPulsation';
% misc.applyLoad = 'planeWave';
misc.model = 'Barrel';
if strcmp(misc.applyLoad,'pointPulsation')
    misc.BC = 'NBC';
else
    misc.BC = 'SHBC';
end
misc.method = {'BEM'};
misc.formulation = {'CCBIE'};
% misc.formulation = {'CCBIE','GBM'};
% misc.formulation = {'GBM'};

varCol = setBarrelParameters(1);
varCol{1}.meshFile = 'createNURBSmesh_Barrel';
f = 1.5e3;             % Frequency
misc.omega = 2*pi*f;

msh.M = 5:6;
msh.M = 1;
msh.degree = 2;
ffp.beta = 0;
msh.parm = [1,2];
msh.parm = 1;
ffp.alpha = (0:0.5:360)*pi/180;

% misc.solveForPtot = true;
misc.solveForPtot = false;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};

prePlot.plot3Dgeometry = 0;
% prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting
err.calculateSurfaceError = true;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;


collectIntoTasks


misc.method = {'KDT'};
misc.solveForPtot = false;
misc.formulation = {'MS1'};

% collectIntoTasks

msh.M = 3;
misc.r_a = 1.25*varCol{1}.R;
pml.t = 0.25*varCol{1}.R;
misc.method = {'PML'};
misc.solveForPtot = false;
misc.formulation = {'GSB'};

% collectIntoTasks


