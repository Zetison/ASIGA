function studies = getTask_Safjan2002tdi_Fig1()
% This study is based on Simpson2014aib and correspond to Figure 6 in Venas2018iao
% Simpson2014aib is available at https://doi.org/10.1016/j.cma.2013.10.026
% Venas2018iao   is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'Safjan2002tdi'; % Simpson sphere
misc.method = {'IENSG'};
% misc.method = {'IE'};

misc.coreMethod = 'IGA';
runTasksInParallel = 0;
misc.progressBars = false;        % Show progress bars for building system matrices
misc.applyLoad = 'Safjan10';
misc.BC = 'NBC';
misc.checkNURBSweightsCompatibility = false;

msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;

c_z = 11;
varCol{1} = struct('media', 'fluid', ...
                   'R',   c_z, ...
                   'c_f', 1500, ...
                   'rho', 1000);
% varCol{1}.c_z = c_z;
% varCol{1}.c_x = sqrt(c_z^2-Upsilon.^2);
% varCol{1}.c_x = varCol{1}.c_z;
msh.meshFile = 'createNURBSmesh_EL';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 1;                 % Wave number for Simpson2014aib
misc.omega = c_f*k;         % Angular frequency

msh.M = 1:7;
msh.M = 7; %8
msh.parm = 1;
ffp.alpha = 0;
ffp.beta = (-90:0.5:90)*pi/180;
ffp.alpha_s = 0;
ffp.beta_s = -pi/2;
prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
prePlot.plotArtificialBndry = false;        % Plot the artificial boundary for the IENSG misc.method
prePlot.plotFullDomain   = 0;        % Plot volumetric domains
prePlot.resolution       = [20,20,0];  % Number of evaluation points in the visualization for each element for each parametric direction
prePlot.view             = [0,0];
prePlot.plotSubsets      = {'xz','Gamma'};
prePlot.plotSubsets      = {'xz'};
prePlot.plotControlPolygon  = 0;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting

misc.computeCondNumber = 0;
err.calculateSurfaceError = 1;
ffp.calculateFarFieldPattern = 1;     % Calculate far field pattern
varCol{1}.refinement = @(M) [0, 2^(M-1)-1, max(2^(M-1)/8-1,0)];
msh.refineThetaOnly = true;
ffp.extraGP = [50,0,0];

msh.degree = 5;
% msh.degree = 2;

warning('off','NURBS:weights')
loopParameters = {'iem.N','iem.p_ie','iem.s_ie','iem.IElocSup', 'iem.IEbasis','misc.method', 'misc.formulation'};

postPlot(1).xname       	= 'iem.N';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType    	= 'loglog';
postPlot(1).lineStyle   	= '*-';
postPlot(1).xLoopName     	= 'iem.N';
postPlot(1).yScale          = 1/100;
% postPlot(1).legendEntries 	= {'misc.method','formulation','M'};
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;

if misc.computeCondNumber
    postPlot(2) = postPlot(1);
    postPlot(2).yname        	= 'cond_number';
    postPlot(2).addCommands   	= [];
end
if ffp.calculateFarFieldPattern
    postPlot(2) = postPlot(1);
    postPlot(2).yname        	= 'TS';
    postPlot(2).xname       	= 'ffp.beta';
    postPlot(2).noXLoopPrms   	= 0;
    postPlot(2).addCommands   	= [];
    postPlot(2).xScale          = 180/pi;
    postPlot(2).yScale          = 1;
    postPlot(2).axisType    	= 'plot';
    postPlot(2).lineStyle   	= '-';
end

iem.IElocSup = false;
iem.s_ie = NaN;
iem.N = 1:19;
iem.N = 1:9;
% iem.N = 15;
iem.p_ie = NaN;
iem.IEbasis	= 'Lagrange';
% iem.IEbasis	= 'Chebyshev';
misc.formulation = {'PGC'};
% misc.formulation = {'PGU'};
iem.boundaryMethod = false;
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
misc.useNeumanProj = 0;
misc.solveForPtot = 0;
iem.N = [1,9];
misc.formulation = {'SL2E'};
postPlot(1).lineStyle = '--';
collectIntoTasks
