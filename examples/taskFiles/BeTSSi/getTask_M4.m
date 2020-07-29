

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4';

method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};

varCol = setM4Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M4';
f = 10e3;             % Frequency
f = 1e3;             % Frequency

M = 3:6;
M = 4;
parm = [2,1];
parm = 2;
degree = 2;
beta = 30*pi/180;
alpha = (0:0.05:360)*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = 0;
prePlot.resolution = [10,10,10];
prePlot.plotNormalVectors = 0;
% prePlot = rmfield(prePlot,'color');
solveForPtot = true;
para.plotResultsInParaview = true;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','formulation','M'};
postPlot(1).subFolderName 	= '../results/M4';
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
    T = readtable('../plotData/refSolutions/M4_HWBC_MS_AS_E30_F1_NTNUFFI.txt','FileType','text', 'HeaderLines',12);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show');
    hold on
end
end
