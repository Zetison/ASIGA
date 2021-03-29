function studies = getTask_Venas2020asi_Figure7(M_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study is based on Venas2020asi (Figure 7, page 316)
% Venas2020asi is available at http://hdl.handle.net/11250/2640443
% Note: To reproduce the data exactly dp_inc = -dp must be used instead of
% the exact expressions for dp_inc computed from e3Dss.m

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

if nargin < 1
    M_0 = 3; % 16
end

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'S1';

misc.method = 'IE';
misc.formulation = 'PGU';
misc.coreMethod = {'SEM','IGA'};
% misc.coreMethod = {'SEM'};
% misc.coreMethod = {'IGA'};
% f = [1e3,1e4];             % Frequency
f = 1e3;             % Frequency
misc.omega = 2*pi*f;

varCol = setS1Parameters('double',1);
msh.meshFile = 'createNURBSmesh_EL';

M = M_0-1; %3:6
iem.N = 3;
misc.N_max = iem.N-1;
ffp.alpha_s = 240*pi/180;
ffp.beta_s = 30*pi/180;
ffp.alpha = (0:0.5:360)*pi/180;
ffp.beta = ffp.beta_s;
err.calculateVolumeError = 1;
err.calculateSurfaceError = 0;
misc.computeCondNumber = true;
misc.clearGlobalMatrices = false;

postPlot(1).xname        	= 'dofsAlg';
postPlot(1).yname        	= 'energyError';
postPlot(1).plotResults  	= 0;
postPlot(1).printResults 	= 0;
postPlot(1).axisType      	= 'loglog';
postPlot(1).lineStyle    	= '*-';
postPlot(1).xScale       	= 1;
postPlot(1).yScale       	= 1;
postPlot(1).xLoopName     	= 'msh.degree';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;

postPlot(2) = postPlot(1);
postPlot(2).axisType = 'loglog';
postPlot(2).xname = 'dofs';
postPlot(2).yname = 'cond_number';

postPlot(3) = postPlot(2);
postPlot(3).yname = 'energyError';
postPlot(3).xname = 'tot_time';

postPlot(4) = postPlot(3);
postPlot(4).xname = 'timeSolveSystem';

postPlot(5) = postPlot(4);
postPlot(5).xname = 'timeBuildSystem';
    

ffp.calculateFarFieldPattern = 0;
prePlot.plot3Dgeometry = 0;
msh.degree = 4:(4+M_0-1);
msh.parm = 2;

loopParameters = {'msh.degree','iem.N','misc.coreMethod','misc.omega'};
collectIntoTasks


