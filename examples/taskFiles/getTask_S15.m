scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'S15';
BC = 'SHBC';
method = {'IE'};
formulation = {'BGU'};

varCol = setS15Parameters();
varCol = varCol(1:3);
varCol{1}.meshFile = 'createNURBSmesh_EL';
% varCol{3}.R_i = 0;
f = 1e3;             % Frequency

M = 1;
parm = 2;
degree = 4;
alpha_s = 240*pi/180;
beta_s = 30*pi/180;
alpha = (0:0.1:360)*pi/180;
beta = 30*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
calculateSurfaceError       = 1;	% Only if scatteringCase == 'Bi'
calculateVolumeError        = 0;	% Only if scatteringCase == 'Bi'
progressBars                = false;        % Show progress bars for building system matrices

prePlot.plot3Dgeometry = 0;
prePlot.resolution = [20,20,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'on';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
% prePlot = rmfield(prePlot,'color');
solveForPtot = false;
prePlot.abortAfterPlotting  = 1;                % Abort simulation after pre plotting

para.plotResultsInParaview = 0;
para.plotDisplacementVectors = false;

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

collectIntoTasks

method = {'BEM'};
solveForPtot = true;
formulation = {'GBM','CCBIE'};
formulation = {'GCBIE'};

% collectIntoTasks
