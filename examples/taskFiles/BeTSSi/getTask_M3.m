function studies = getTask_M3()
% This study is based on Figure B.14 in Venas2019asi
% Venas2019asi is available at http://hdl.handle.net/11250/2640443

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.applyLoad = 'planeWave'; % Set load. I.e.: 'planeWave', 'radialPulsation', 'pointPulsation', 'SimpsonTorus'
% misc.applyLoad = 'pointPulsation';

misc.scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M3';
misc.BC = 'SHBC';
% misc.BC = 'NBC';
misc.method = {'IENSG'};
misc.formulation = 'BGU';
misc.checkNURBSweightsCompatibility = false;
err.calculateSurfaceError = 0;

prePlot.plot2Dgeometry = 0;
prePlot.plot3Dgeometry = 0;
prePlot.resolution = [80,60,0];
prePlot.elementBasedSamples = 0;
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotGeometryInfo = 1;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting
prePlot.pngResolution = '-r800';
% prePlot.LineWidth           = 0.5;         % Width of lines

varCol = setM3Parameters(1);
varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-1)/8-1, 2^(M-1)/4-1];
msh.meshFile = 'createNURBSmesh_M3';
msh.degree = 2;
% msh.degree = 5; % use and odd degree > 4 if parm = 2 and misc.method = 'IENSG' (singular evaluation at poles not yet implemented for 'IENSG')
msh.parm = 1;
msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
msh.explodeNURBS = ~prePlot.plot2Dgeometry;
msh.M = 2;
f = 1e3;             % Frequency
misc.omega = 2*pi*f;

ffp.beta = 0;
ffp.alpha = (0:0.1:360)*pi/180;

warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation'};
misc.solveForPtot = false;

para.plotResultsInParaview = false;

postPlot(1).xname       	= 'alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_(i_study);

% postPlot(2)                 = postPlot(1);
% postPlot(2).xname       	= 'nepw';
% postPlot(2).yname        	= 'surfaceError';
% postPlot(2).plotResults  	= 0;
% postPlot(2).printResults 	= 0;
% postPlot(2).axisType        = 'loglog';
% postPlot(2).lineStyle   	= '-*';
% postPlot(2).xLoopName     	= 'msh.M';
% postPlot(2).fileDataHeaderX	= [];
% postPlot(2).noXLoopPrms   	= 1;
% postPlot(2).xScale          = 1;
% postPlot(2).addCommands   	= [];

if strcmp(misc.scatteringCase,'MS')
    para.i_MS = find(abs(ffp.alpha - 240*pi/180) < 20*eps);
    postPlot(1).xlim            = [0,180];
    postPlot(1).ylim            = [-60,40];
else
    postPlot(1).xlim            = [0,360];
    postPlot(1).ylim            = [-40,50];
    ffp.beta_s = 0;
    ffp.alpha_s = 240*pi/180;
end

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
misc.formulation = {'MS1'};
% misc.coreMethod = 'linear_FEM';
msh.M = 5; % 5

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEM simulation
misc.solveForPtot = strcmp(misc.applyLoad,'planeWave');
misc.method = {'BEM'};
misc.formulation = {'CCBIE','CBM'};
misc.formulation = {'CCBIE'};
% misc.formulation = {'GBM'};
msh.M = 5:6;
msh.degree = 4:5;
msh.parm = 2;
% msh.M = 1;
loopParameters = {'misc.formulation','msh.M','msh.degree','misc.method','misc.omega'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
misc.solveForPtot = false;
misc.formulation = '';
misc.progressBars = 0;
err.err.calculateSurfaceError = 0;
misc.computeCondNumber = false;
ffp.plotFarField = 1;
msh.M = 1;
misc.applyLoad = 'planeWave';
msh.degree = 2;
msh.parm = 1;
rt.N = 6; % 6
msh.x_0 = -[varCol{1}.L/2+(varCol{1}.R2-varCol{1}.R1)/2,0,0]; 
% rt.N = 1; % 6
warning('off','RT:limitations')
loopParameters = {'msh.parm','msh.M','misc.method','rt.N'};
% collectIntoTasks

msh.x_0 = [0,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PML simulation
misc.method = {'PML'};
misc.formulation = {'GSB'};
msh.M = 4:6;
% msh.M = 5;
msh.parm = 1;

pml.t = 0.25*varCol{1}.R2;         % thickness of PML
misc.r_a = 1.25*varCol{1}.R2;         % thickness of PML

para.plotResultsInParaview	 = 0;	% Only if misc.scatteringCase == 'Bi'
para.extraXiPts              = 'round(60/2^(M-1))';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(60/2^(M-1))';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(60/2^(M-1))';  % Extra visualization points in the zeta-direction per element
para.plotTimeOscillation     = 0;
para.plotFullDomain          = 0;
para.plotSubsets             = {'Gamma','Gamma_a','xy','xz'}; % Plot subsets (i.e. the artificial boundary Gamma_a) in paraview

loopParameters = {'misc.formulation','msh.degree','msh.M','misc.method','misc.omega'};
% collectIntoTasks


function addCommands_(i_study)
if i_study == 1
    if 1
        T = readtable('miscellaneous/refSolutions/M3_BEM_IGA_CBM_M7_f1000_N10_TSVSalpha.txt', ...
                        'FileType','text','CommentStyle','%');
        x = T.alpha;
        y = T.TS;
        plot(x,y,'DisplayName','Reference solution f = 1000Hz')
    else
        T = readtable('miscellaneous/refSolutions/M3_SHBC_BI_alpha240.txt', ...
                        'FileType','text','headerLines',6);
        x = T.Var1;
        y = T.Var2;
        plot(x,y,'DisplayName','Reference solution f = 1000Hz')
        
        T = readtable('miscellaneous/refSolutions/M3_HWBC_BI_A240_E0_F1_FOI.txt', ...
                        'FileType','text','headerLines',6);
        x = T.Var1;
        y = T.Var2;
        plot(x,y,'DisplayName','FOI')
    end
    legend('off');
    legend('show','Interpreter','latex');
end


