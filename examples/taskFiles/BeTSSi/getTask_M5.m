function studies = getTask_M5()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = {'M5A', 'M5B'}; % BeTSSi misc.model 5A and BeTSSi misc.model 5B
% misc.model = {'M5A'}; % BeTSSi misc.model 5A and BeTSSi misc.model 5B
misc.BC = 'SHBC';
misc.method = {'BEM'};
misc.formulation = {'CBM','CCBIE'};
misc.formulation = {'CCBIE'};

varCol = setM5Parameters();
msh.meshFile = 'createNURBSmesh_M5';
f = 1e3;             % Frequency
misc.omega = 2*pi*f;
subFolderName = 'M5';

msh.M = 4:5;
msh.M = 3;
msh.parm = [1,2];

ffp.degree = 2;
ffp.beta = 0;
ffp.alpha = (0:0.5:360)*pi/180;
ffp.alpha_s = 60*pi/180;
ffp.beta_s = 0*pi/180;
misc.solveForPtot = true;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation','misc.model'};

prePlot.plot3Dgeometry = 0;
prePlot.view = [7,30];
% prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 1;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting


postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.model', 'misc.method','msh.parm','misc.formulation','msh.M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

collectIntoTasks


misc.method = {'KDT'};
misc.solveForPtot = false;
misc.formulation = {'MS1'};

collectIntoTasks



function addCommands_(i_study)
T = readtable('miscellaneous/refSolutions/M5A_HWBC_BI_A60_E0_F1.txt','FileType','text', 'HeaderLines',7);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','M5A: Reference solution f = 1000Hz')
T = readtable('miscellaneous/refSolutions/M5B_HWBC_BI_A60_E0_F1.txt','FileType','text', 'HeaderLines',7);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','M5B: Reference solution f = 1000Hz')
legend('off');
legend('show');
hold on