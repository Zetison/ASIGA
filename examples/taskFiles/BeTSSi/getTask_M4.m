function studies = getTask_M4(M_0)
% This study is based on Figure 4 page 272 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

if nargin < 1
    M_0 = 6; % 
end

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M4';

misc.method = {'BEM'};
misc.formulation = {'CBM','CCBIE'};
misc.formulation = {'CCBIE'};
misc.formulation = {'CCBIE','GBM'};

varCol = setM4Parameters(1);
msh.meshFile = 'createNURBSmesh_M4';
% f = 3e3;             % Frequency
f = 10e3;             % Frequency
misc.omega = 2*pi*f;
msh.M = 5:M_0;
% msh.M = 1;
msh.parm = [2,1];
msh.parm = 2;
msh.degree = 2;
ffp.beta = 30*pi/180;
ffp.alpha = (0:0.05:360)*pi/180;
misc.solveForPtot = true;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};
prePlot.plot3Dgeometry = 0;
prePlot.view = [120,20];
% prePlot.resolution = [20,20,20];
prePlot.plotNormalVectors = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting


para.plotResultsInParaview = false;

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= 0;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
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
if i_study == 1
    f = 10000;
    if f == 1000
        T = readtable('miscellaneous/refSolutions/M4_HWBC_MS_AS_E30_F1.txt','FileType','text', 'HeaderLines',7);
        x = T.Var1*pi/180;
        y = T.Var2;
    else
        T = readtable('miscellaneous/refSolutions/M4_M7_parm1_f10000_BEM_CCBIE_TSVSalpha.txt','FileType','text','CommentStyle','%');
        x = T.alpha*pi/180;
        y = T.TS;
    end
    rlim([-70,40])
    polarplot(x,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show');
    hold on
end
