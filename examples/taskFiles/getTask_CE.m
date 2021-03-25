function studies = getTask_CE()
% This study is based on Figure 6 and 7 in Makarov1998ais
% Makarov1998ais is available at https://doi.org/10.1121/1.421238

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'CE';
misc.method = {'BEM'};
misc.formulation = {'CCBIE'};
misc.formulation = {'CCBIE'};
% misc.formulation = {'GBM'};

varCol{1} = struct('media', 'fluid', ...
                   'R', 1, ...
                   'c_f', 1500, ...
                   'rho', 1000);
msh.meshFile = 'createNURBSmesh_CE';
misc.omega = 62.8*varCol{1}.c_f;             % angular frequency

misc.omega = 10000;             % angular frequency

misc.r_a = 1.25*varCol{1}.R;         % thickness of PML
msh.M = 6:7;
msh.M = 4;
msh.parm = 1;
msh.degree = 2;
ffp.beta = 0;
ffp.alpha = (0:0.05:360)*pi/180;
ffp.beta_s = 0;
ffp.alpha_s = pi/4;
misc.progressBars = false;        % Show progress bars for building system matrices
misc.solveForPtot = true;
warning('off','NURBS:weights')

loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.resolution = [10,10,0];
prePlot.view = [135,30];
prePlot.plotParmDir = 0;
prePlot.plotControlPolygon  = 0;                 % Plot the control polygon for the NURBS mesh
prePlot.plotNormalVectors = 0;
prePlot.abortAfterPlotting  = 1;                % Abort simulation after pre plotting
prePlot.plotGeometryInfo    = 1;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

% collectIntoTasks

misc.method = {'KDT'};
misc.solveForPtot = false;
misc.formulation = {'MS1'};

collectIntoTasks

misc.method = {'IE'};
misc.formulation = {'BGU'};

collectIntoTasks

pml.t = 0.25*varCol{1}.R;         % thickness of PML
varCol{1}.refinement = @(M) [0, 0, 0, 2^(M-1)/2-1];
misc.method = {'PML'};
misc.formulation = {'GSB'};

collectIntoTasks

msh.M = 5; % 7
para.plotResultsInParaview = 1;
para.extraXiPts              = 'round(20/2^(M-1))';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(20/2^(M-1))';  % Extra visualization points in the eta-direction per element

% collectIntoTasks


