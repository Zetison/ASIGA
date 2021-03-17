function studies = getTask_M3()
% This study is based on Figure B.14 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M3';
BC = 'SHBC';
misc.method = {'IENSG'};
formulation = 'BGU';

varCol = setM3Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M3';
f = 1e3;             % Frequency

M = 4:5;
% M = 1;
degree = 2;
% degree = 5; % use and odd degree > 4 if parm = 2 and misc.method = 'IENSG' (singular evaluation at poles not yet implemented for 'IENSG')
parm = 1;
beta = 0;
alpha = (0:0.1:360)*pi/180;

warning('off','NURBS:weights')
loopParameters = {'M','parm','f','misc.method','formulation'};
prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 0;                % Abort simulation after pre plotting
solveForPtot = false;

para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'misc.method','parm','f','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).xlim            = [0,180];
postPlot(1).ylim            = [-60,40];
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
formulation = {'MS1'};

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
solveForPtot = true;
misc.method = {'BEM'};
formulation = {'CCBIE','CBM'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};
BC = 'SHBC';
loopParameters = {'formulation','M','misc.method','f'};
collectIntoTasks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
solveForPtot = false;
formulation = '';
err.err.calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
M = 2;
misc.applyLoad = 'planeWave';
N = 4:6;
r = 10;

loopParameters = {'parm','M','misc.method','N'};
collectIntoTasks


function addCommands_(i_study)
if i_study == 1
    T = readtable('miscellaneous/refSolutions/M3_BEM_IGA_CBM_M7_f1000_N10_TSVSalpha.txt', ...
                    'FileType','text', 'HeaderLines',1,'CommentStyle','%');
    x = T.Var1;
    y = T.Var2;
    polarplot(x*pi/180,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show','Interpreter','latex');
end


