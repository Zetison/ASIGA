%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study is based on Makarov1998ais and correspond to Figure 6 and 7 in Venas2018iao
% Makarov1998ais is available at https://doi.org/10.1121/1.421238

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'CE';

method = {'BEM'};
formulation = {'CCBIE'};

varCol{1} = struct('media', 'fluid', ...
                   'R', 1, ...
                   'c_f', 1500, ...
                   'rho', 1000);
varCol{1}.meshFile = 'createNURBSmesh_CE';
f = 62.8*varCol{1}.c_f/(2*pi);             % Frequency

M = 5:6;
parm = 2;
degree = 4;
beta = 0;
alpha = (0:0.05:360)*pi/180;
beta_s = 0;
alpha_s = pi/4;

quadMethodBEM = 'Simpson';
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = 1;
prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 1;
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'M';
postPlot(1).legendEntries 	= {'method','parm','f','formulation','M','alpha_s'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

collectIntoTasks

M = 6;
para.plotResultsInParaview = 1;
para.plotDisplacementVectors = false;

collectIntoTasks

method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks