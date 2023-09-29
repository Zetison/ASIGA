function studies = getTask_Venas2019asi_FigureA2A3(M_0)
% This study is based on Figure A.2 and Figure A.3 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

if nargin < 1
    M_0 = 4; 
end

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M1';

misc.method = {'BEM'};
misc.formulation = {'CBM','CCBIE'};
misc.formulation = {'CCBIE'};
% misc.formulation = {'CCBIE','GBM'};

varCol = setM1Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M1';
f = [1e2,1e3];             % Frequency
f = 1e2;             % Frequency
misc.omega = 2*pi*f;

msh.M = 4:5;
msh.M = M_0;
msh.degree = 2;
ffp.beta = 0;
ffp.alpha = (0:0.5:360)*pi/180;
misc.solveForPtot = true;

loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};
prePlot.plot3Dgeometry = 0;
% prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.view = [198,10];
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting

para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','msh.parm','misc.omega','misc.formulation','msh.M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% collectIntoTasks

misc.solveForPtot = false;
misc.method = {'IENSG'};
misc.formulation = {'BGU'};
iem.N = [3,5];
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation','iem.N'};
collectIntoTasks

misc.method = {'KDT'};
misc.formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
for f = 1e2 %[1e2, 1e3]
    T = readtable(['miscellaneous/refSolutions/M1_BEM_IGA_GBM_M6_NNaN_f' num2str(f) '_TSVSalpha.txt'], ...
                    'FileType','text','CommentStyle','%');
    x = T.alpha;
    y = T.TS;
    polarplot(x*pi/180,y,'DisplayName',['Reference solution f = ' num2str(f) 'Hz'])
    legend('off');
    legend('show','Interpreter','latex');
end