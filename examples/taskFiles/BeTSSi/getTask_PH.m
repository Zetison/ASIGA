function studies = getTask_PH()
% This study is based on Schneider2003asb Figure 3
% Schneider2003asb is available at https://www.researchgate.net/publication/27256874

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'MS';
misc.model = 'PH'; % Simpson sphere
misc.BC = 'SHBC';
misc.method = {'IENSG','IE'};
misc.method = {'IE'};
misc.formulation = 'BGU';
misc.coreMethod = 'IGA';
misc.computeCondNumber = 0;

varCol = setPHParameters(1);
msh.meshFile = 'createNURBSmesh_PH';
f = 2e2; %[2e2 1e3 4e3];             % Frequency
misc.omega = 2*pi*f;

msh.M = 3:4;
% msh.M = 1;

% ffp.alpha = (0:0.1:360)*pi/180;
ffp.alpha = (0:1:360)*pi/180;

prePlot.plot3Dgeometry = 0;
prePlot.plotControlPolygon = false;
prePlot.plotNormalVectors = false;
% prePlot.resolution = [20,20,0];
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting

msh.degree = 2;

iem.N = 3;

loopParameters = {'msh.M','misc.formulation','misc.method'};
msh.parm = 1;

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).legendEntries 	= {'misc.method','misc.formulation','msh.M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
misc.formulation = 'MS1';
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
% formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
misc.formulation = {'CCBIE'};
misc.formulation = {'CCBIE','GBM'};
misc.solveForPtot = true;
% collectIntoTasks


function addCommands_(i_study)
T = readtable('miscellaneous/refSolutions/PH_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',12);
x = T.Var1;
y = T.Var2;
plot(x,y,'DisplayName','Reference solution')
legend('off');
legend('show');
hold on