scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M1';

method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};

varCol = setM1Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M1';
f = 1e3;             % Frequency

M = 3;
degree = 2;
beta = 0;
alpha = (0:0.05:360)*pi/180;

alpha_s = [240, 300]*pi/180;
beta_s = 0;

loopParameters = {'M','parm','f','method','formulation','alpha_s'};
prePlot.plot3Dgeometry = 0;
prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 1;
prePlot.axis = 'on';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
% prePlot = rmfield(prePlot,'color');
solveForPtot = true;

para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','f','formulation','M','alpha_s'};
postPlot(1).subFolderName 	= '../results/M1';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% collectIntoTasks


M = 2;
method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
if i_study == 1
    T = readtable('../plotData/refSolutions/M1_HWBC_BI_0_1.txt','FileType','text', 'HeaderLines',7);
    x = T.Var1;
    y240 = T.Var2;
    y300 = T.Var3;
    plot(x,y240,'DisplayName','Reference solution f = 1000Hz, at $$\alpha_s = 240^\circ$$')
    hold on
    plot(x,y300,'DisplayName','Reference solution f = 1000Hz, at $$\alpha_s = 300^\circ$$')
    legend('off');
    legend('show','Interpreter','latex');
end
end