function studies = getTask_PH()
% This study is based on Schneider2003asb Figure 3
% Schneider2003asb is available at https://www.researchgate.net/publication/27256874

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'MS';
misc.model = 'PH'; % Simpson sphere
BC = 'SHBC';
misc.method = {'IENSG','IE'};
% misc.method = {'IE'};
formulation = 'BGU';
misc.coreMethod = 'IGA';
computeCondNumber = 0;
runTasksInParallel = 0;

varCol = setPHParameters(1);
% varCol{1}.meshFile = 'createNURBSmesh_PH';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
f = 2e2; %[2e2 1e3 4e3];             % Frequency

M = 3:4;
% M = 1;

% alpha = (0:0.05:360)*pi/180;
alpha = (0:1:360)*pi/180;

prePlot.plot3Dgeometry = 1;
prePlot.plotControlPolygon = false;
prePlot.plotNormalVectors = false;
% prePlot.resolution = [20,20,0];
prePlot.abortAfterPlotting = false;                % Abort simulation after pre plotting

degree = 2;

N = 3;

loopParameters = {'M','formulation','misc.method'};
parm = 1;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
formulation = 'MS1';
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
% formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};
solveForPtot = true;
collectIntoTasks


function addCommands_(i_study)
if i_study == 1
    T = readtable('miscellaneous/refSolutions/PH_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',12);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName','Reference solution')
    legend('off');
    legend('show');
    hold on
end