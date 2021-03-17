function studies = getTask_M1()
% This study is based on Figure A.2 and Figure A.3 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M1';

misc.method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};

varCol = setM1Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M1';
f = [1e2,1e3];             % Frequency
% f = 1e2;             % Frequency

M = 4:5;
% M = 1;
degree = 2;
beta = 0;
alpha = (0:0.5:360)*pi/180;
solveForPtot = true;

loopParameters = {'M','parm','f','misc.method','formulation','alpha_s'};
prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.view = [198,10];
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 0;                % Abort simulation after pre plotting

para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','parm','f','formulation','M','alpha_s'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

collectIntoTasks

solveForPtot = false;
misc.method = {'IENSG'};
formulation = {'BGU'};
N = 3;
collectIntoTasks

misc.method = {'KDT'};
formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
if i_study == 1
    for f = 1e2 %[1e2, 1e3]
        T = readtable(['miscellaneous/refSolutions/M1_BEM_IGA_GBM_M6_NNaN_f' num2str(f) '_TSVSalpha.txt'], ...
                        'FileType','text', 'HeaderLines',1,'CommentStyle','%');
        x = T.Var1;
        y = T.Var2;
        polarplot(x*pi/180,y,'DisplayName',['Reference solution f = ' num2str(f) 'Hz'])
        legend('off');
        legend('show','Interpreter','latex');
    end
end