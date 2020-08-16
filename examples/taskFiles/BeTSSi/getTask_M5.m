scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

% model = {'M5A', 'M5B'}; % BeTSSi model 5A and BeTSSi model 5B
model = {'M5A'}; % BeTSSi model 5A and BeTSSi model 5B
BC = 'SHBC';
method = {'BEM'};
formulation = {'CBM','CCBIE'};
formulation = {'CCBIE'};

varCol = setM5Parameters();
varCol{1}.meshFile = 'createNURBSmesh_M5';
varCol{1}.type = 'B'; % A or B (M5A or M5B)
f = 1e3;             % Frequency
subFolderName = 'M5';

M = 1;
degree = 2;
beta = 0;
parm = 2;
alpha = (0:0.5:360)*pi/180;
alpha_s = 60*pi/180;
beta_s = 0*pi/180;
quadMethodBEM = 'Simpson';
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};

prePlot.plot3Dgeometry = 1;
prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;


postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'model', 'method','parm','formulation','M'};
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
    T = readtable('miscellaneous/refSolutions/M5B_HWBC_BI_A60_E0_F1.txt','FileType','text', 'HeaderLines',7);
    x = T.Var1;
    y = T.Var2;
    plot(x,y,'DisplayName','Reference solution f = 1000Hz')
    legend('off');
    legend('show');
    hold on
end
end
