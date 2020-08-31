function studies = getTask_M4()
% This study is based on Figure 4 page 273 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4';

method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};

varCol = setM4Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M4';
% f = 3e3;             % Frequency
f = 1e3;             % Frequency

M = 4:5;
% M = 1;
parm = [2,1];
parm = 2;
degree = 2;
beta = 30*pi/180;
alpha = (0:0.05:360)*pi/180;
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = 1;
prePlot.view = [120,20];
% prePlot.resolution = [20,20,20];
prePlot.plotNormalVectors = 0;
prePlot.abortAfterPlotting = 0;                % Abort simulation after pre plotting


para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

collectIntoTasks


method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
if i_study == 1
    T = readtable('miscellaneous/refSolutions/M4_HWBC_MS_AS_E30_F1.txt','FileType','text', 'HeaderLines',7);
    x = T.Var1*pi/180;
    y = T.Var2;
    polarplot(x,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show');
    hold on
end
