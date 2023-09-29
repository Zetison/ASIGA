function studies = getTask_M2()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M2';

misc.method = {'BEM'};
misc.formulation = {'CCBIE'};
misc.formulation = {'CCBIE','GBM'};

varCol = setM2Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M2';
f = 1e3;
misc.omega = 2*pi*f;             % Frequency

msh.M = 4:5;
% msh.M = 5;
msh.degree = 2;
ffp.eta = 0;
ffp.alpha = (0:0.1:360)*pi/180;
misc.solveForPtot = true;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};
prePlot.plot3Dgeometry = 0;
% prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).legendEntries 	= {'misc.method','msh.parm','misc.formulation','msh.M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% collectIntoTasks

misc.method = {'KDT'};
misc.solveForPtot = false;
misc.formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
T = readtable('../../OneDrive/work/matlab/plotData/refSolutions/M2_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',7);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','Reference solution f = 1000Hz')
legend('off');
legend('show');
hold on
% 
