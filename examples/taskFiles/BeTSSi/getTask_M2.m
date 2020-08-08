scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M2';

method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};

varCol = setM2Parameters(1);
varCol{1}.meshFile = 'createNURBSmesh_M2';
f = 1e3;             % Frequency

M = 3:4;
degree = 2;
beta = 0;
alpha = (0:0.5:360)*pi/180;

pctRunOnAll warning('off', 'BEM:recursion')

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = 0;
prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'on';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
% prePlot = rmfield(prePlot,'color');
solveForPtot = true;
para.plotResultsInParaview = true;
prePlot.abortAfterPlotting  = 1;                % Abort simulation after pre plotting

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% collectIntoTasks


M = 4;
method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks

function addCommands_(i_study)
if i_study == 1
    T = readtable('../../OneDrive/work/matlab/plotData/refSolutions/M2_HWBC_MS_AS_E0_F1.txt','FileType','text', 'HeaderLines',7);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show');
    hold on
end
end


