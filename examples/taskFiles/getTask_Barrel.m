function studies = getTask_Barrel()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

% misc.applyLoad = 'pointPulsation';
misc.applyLoad = 'planeWave';
misc.model = 'Barrel';
if strcmp(misc.applyLoad,'pointPulsation')
    misc.BC = 'NBC';
else
    misc.BC = 'SHBC';
end
misc.method = {'BEM'};
misc.formulation = {'CCBIEC'};
% misc.formulation = {'CCBIE','GBM'};
% misc.formulation = {'GBM'};

varCol = setBarrelParameters(1);
R = varCol{1}.R;
L = varCol{1}.L;
misc.r_a = 1.25*R;
pml.t = 0.25*R;
pml.refinement = @(M) round((2^(M-1)-1)*pml.t/(R*2*pi/3));
msh.meshFile = 'createNURBSmesh_Barrel';
% varCol{1}.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
varCol{1}.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
k = 100;
% k = 10;
misc.omega = k*varCol{1}.c_f;
msh.refineThetaOnly = true;
msh.pmlFill = true;
msh.M = 8:9;
msh.M = 1; % 8
msh.degree = 2;
msh.parm = [1,2];
msh.parm = 1;
ffp.alpha_s = pi;
ffp.beta_s = 0;
ffp.extraGP = [7,0,0];

ffp.beta = 0;
ffp.alpha = (0:0.05:360)*pi/180;

% misc.solveForPtot = true;
misc.solveForPtot = false;
misc.checkNURBSweightsCompatibility = false; 
warning('off','NURBS:weights')
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation','misc.applyLoad','msh.pmlFill'};

prePlot.plot3Dgeometry = 0;
prePlot.resolution = [100,0,20];
% prePlot.resolution = [0,0,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting
prePlot.coarseLinearSampling = prePlot.plotParmDir;
prePlot.plotSubsets          = {'xz','Gamma'};
prePlot.plotSubsets          = {'Gamma'};
prePlot.plotFullDomain       = 0;
prePlot.view = [0,0];
prePlot.camproj = 'orthographic';
prePlot.useCamlight = false;

err.calculateSurfaceError = strcmp(misc.applyLoad,'pointPulsation');

postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'TS';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'polar';
postPlot(1).lineStyle   	= '-';
postPlot(1).xLoopName     	= 'msh.M';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 180/pi;

para.plotFullDomain          = false;
para.plotResultsInParaview	 = 1;
para.extraXiPts              = 'max(round(2^(M-6)-1),1)';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'max(round(2^(M-6)-1),1)';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'max(round(2^(M-6)-1),1)';  % Extra visualization points in the zeta-direction per element
para.plotSubsets             = {'xz','Gamma','Gamma_a'}; % Plot (surface) subsets (i.e. the artificial boundary Gamma_a) in paraview 
                                                    % (examples include: 'Gamma','Gamma_a','yz','xz','xy','innerCoupling','outerCoupling','outer','inner','homDirichlet')
msh.pmlFill = false; % Use rounded corners for PML domain
misc.method = {'PML'};
misc.solveForPtot = false;
misc.formulation = {'GSB'};

% collectIntoTasks

iem.boundaryMethod = 0;
misc.method = {'IENSG'};
misc.formulation = {'BGC'};
% collectIntoTasks

para.plotResultsInParaview = 0;
msh.M = 6:7;
% msh.M = 1;
msh.parm = 2;
msh.degree = 2;
msh.refineThetaOnly = false;
misc.method = {'BEM'};
misc.formulation = {'CCBIE'};
for applyLoad = {'pointPulsation','planeWave'}
    err.calculateSurfaceError = strcmp(applyLoad,'pointPulsation');
    misc.applyLoad = applyLoad{1};
    if strcmp(misc.applyLoad,'pointPulsation')
        misc.BC = 'NBC';
    else
        misc.BC = 'SHBC';
    end
    misc.solveForPtot = strcmp(misc.applyLoad,'planeWave');
    collectIntoTasks
end

%% Do frequency analysis for BEM
postPlot(1).xname       	= 'varCol{1}.kL';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= false;
postPlot(1).axisType        = 'semilogy';
postPlot(1).lineStyle   	= '-';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale          = 1;

msh.degree = 2;
misc.applyLoad = 'pointPulsation';
ffp.alpha = 0;
misc.method = {'BEM'};
misc.formulation = {'CCBIE','CHBIE','CBM','CCBIEC'};
misc.formulation = {'CCBIE'};
misc.applyLoad = 'pointPulsation';
err.calculateSurfaceError = strcmp(misc.applyLoad,'pointPulsation');
misc.solveForPtot = ~strcmp(misc.applyLoad,'pointPulsation');
misc.BC = 'NBC';
msh.M = 4;

%% Run convergence analysis
misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
postPlot(1).xname       	= 'surfDofs';
postPlot(1).noXLoopPrms   	= 1;
msh.degree = 2;
k = 10;
msh.parm = 1;
misc.solveForPtot = false;
msh.refineThetaOnly = true;
misc.omega = k*varCol{1}.c_f;
misc.method = {'BA'};
misc.solveForPtot = false;
misc.formulation = {'SL2E'};
msh.M = 1:8; % 8
para.plotResultsInParaview	 = 0;
loopParameters = {'msh.M','msh.parm','misc.omega','misc.method','misc.formulation','msh.pmlFill'};

% collectIntoTasks

msh.pmlFill = [0,1];
% msh.pmlFill = 1;
misc.method = {'PML'};
misc.solveForPtot = false;
misc.formulation = {'GSB'};

% collectIntoTasks

