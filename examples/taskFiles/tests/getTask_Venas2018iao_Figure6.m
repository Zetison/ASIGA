function studies = getTask_Venas2018iao_Figure6(M_0)
% This study is based on Simpson2014aib and correspond to Figure 6 in Venas2018iao
% Simpson2014aib is available at https://doi.org/10.1016/j.cma.2013.10.026
% Venas2018iao   is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)

if nargin < 1
    M_0 = 3;
    plotResults = true;
else
    plotResults = false;
end
counter = 1;
studies = cell(0,1);
getDefaultTaskValues

%% IE simulation
misc.scatteringCase = 'BI';
misc.model = 'SS'; % Simpson sphere
misc.method = {'IE'};
misc.BC = 'SHBC';
misc.formulation = {'PGU','PGC','BGU','BGC'};
% misc.formulation = {'BGU'};
misc.coreMethod = 'IGA';
misc.computeCondNumber = 0;
runTasksInParallel = 0;
misc.progressBars = false;        % Show progress bars for building system matrices

warning('off','NURBS:weights')

varCol = setSSParameters(1);
msh.meshFile = 'createNURBSmesh_EL';
c_f = varCol{1}.c_f;   % Speed of sound in outer fluid
k = 2;                 % Wave number for Simpson2014aib
misc.omega = c_f*k;         % Angular frequency

msh.M = M_0; % 3:4
msh.degree = 3;
msh.parm = 1;
msh.autoRefine = true;

ffp.plotFarField = false; 
ffp.alpha_s = pi;
ffp.beta_s = 0;  
ffp.alpha = (0:0.5:360)*pi/180;
ffp.r = 5; % radii for near-field evaluation

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
prePlot.abortAfterPlotting  = true;       % Abort simulation after pre plotting
% prePlot.resolution = [20,20,0];

err.calculateSurfaceError = 1;

iem.N = 4;

loopParameters = {'msh.M','misc.method','misc.formulation'};

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'error_pAbs';
postPlot(1).plotResults  	= plotResults;
postPlot(1).printResults 	= false;
postPlot(1).axisType    	= 'semilogy';
postPlot(1).lineStyle   	= '-';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms     = 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).yScale          = 1/100;
postPlot(1).addCommands   	= [];
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_();

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BA simulation
misc.method = {'BA'};
misc.solveForPtot = 1;
% misc.formulation = {'SL2E','VL2E'};
misc.formulation = {'SL2E'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.method = {'BEM'};
msh.M = 1:M_0; % 1:3
% misc.formulation = {'CCBIE','CHBIE','CBM','GCBIE','GHBIE','GBM'};
misc.formulation = {'CCBIE'};
% misc.formulation = {'GCBIE','GHBIE','GBM'};
bem.colBEM_C0 = 0;
misc.solveForPtot = true;
collectIntoTasks

function addCommands_()
    
error_simpson = importdata('miscellaneous/refSolutions/Fig17_M1.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson, M=1');
error_simpson = importdata('miscellaneous/refSolutions/Fig17_M2.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson, M=2');
error_simpson = importdata('miscellaneous/refSolutions/Fig17_M3.csv');
loglog(180/pi*error_simpson(:,1),error_simpson(:,2),'*','DisplayName','Simpson, M=3');
legend('off');
legend('show','Interpreter','latex');