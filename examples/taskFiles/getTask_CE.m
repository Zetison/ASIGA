function studies = getTask_CE()
% This study is based on Makarov1998ais and correspond to Figure 6 and 7 in Venas2018iao
% Makarov1998ais is available at https://doi.org/10.1121/1.421238

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'CE';
method = {'BEM'};
formulation = {'CCBIE'};
formulation = {'CCBIE','GBM'};

varCol{1} = struct('media', 'fluid', ...
                   'R', 1, ...
                   'c_f', 1500, ...
                   'rho', 1000);
varCol{1}.meshFile = 'createNURBSmesh_CE';
f = 62.8*varCol{1}.c_f/(2*pi);             % Frequency

M = 6:7;
% M = 6;
parm = 2;
degree = 4;
beta = 0;
alpha = (0:0.05:360)*pi/180;
beta_s = 0;
alpha_s = pi/4;
progressBars = false;        % Show progress bars for building system matrices
solveForPtot = true;

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [20,20,0];
prePlot.view = [135,30];
prePlot.elementBasedSamples = 0;
prePlot.plotParmDir = 0;
prePlot.plotControlPolygon  = 0;                 % Plot the control polygon for the NURBS mesh
prePlot.plotNormalVectors   = 0;                % Plot the normal vectors for the NURBS mesh
prePlot.abortAfterPlotting  = 0;                % Abort simulation after pre plotting

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

method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};

collectIntoTasks

method = {'BEM'};
solveForPtot = true;
formulation = {'GBM'};
M = 7;
para.plotResultsInParaview = 1;
para.plotDisplacementVectors = false;
para.extraXiPts              = 'round(40/2^(M-1))';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(40/2^(M-1))';  % Extra visualization points in the eta-direction per element

collectIntoTasks