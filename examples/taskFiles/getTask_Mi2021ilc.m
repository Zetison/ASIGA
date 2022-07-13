function studies = getTask_Mi2021ilc()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 19 to 22 in Mi2021ilc
% Mi2021ilc is available at https://doi.org/10.1016/j.cma.2021.113925

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'Mi2021ilc';
% noCoresToUse = 4;
msh.explodeNURBS = 1;   % Create patches from all C^0 interfaces

msh.meshFile = 'createNURBSmesh_EL';
msh.nonLinearParam = true;
msh.parm = 1;
misc.checkNURBSweightsCompatibility = 0;
prePlot.plotGeometryInfo    = 0;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
err.calculateVolumeError    = true;
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}

prePlot.plotFullDomain   = 1;        % Plot volumetric domains
prePlot.view             = [0,90];
prePlot.plotSubsets      = {'xy'};
% prePlot.plotSubsets      = {};
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.plotControlPolygon = 1;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting = true;       % Abort simulation after pre plotting
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
warning('off','NURBS:weights')

% postPlot(1).xname       	= 'nepw';
postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'abs_p';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 1;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale       	= 180/pi;
postPlot(1).legendEntries 	= {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','varCol{1}.k'};
postPlot(1).addCommands   	= @(study,i_study,studies) addCommands_();
postPlot(2) = postPlot(1);
postPlot(2).yname        	= 'error_pAbs';
postPlot(2).axisType        = 'semilogy';
postPlot(2).addCommands   	= @(study,i_study,studies) addCommands2_();

msh.meshFile = 'createNURBSmesh_EL';
% msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
msh.refineThetaOnly = false;

connectedParameters = {{'misc.method','misc.formulation'},{'pml.sigmaType','pml.n'}};
% misc.method = {'PML','IE'};
% misc.formulation = {'GSB','BGU'};
misc.method = {'PML'};
misc.formulation = {'GSB'};
misc.coreMethod = {'IGA'};

R = 1;
r = 0.2;
c_f = 340;
msh.M = 2:4; % 4
% msh.M = 1; % 4
varCol{1}.media = 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
varCol{1}.R = r;
varCol{1}.c_f = c_f;
varCol{1}.rho = 1.225;

kr = [pi,2*pi,3*pi];
% kr = 3*pi;
% kr = pi;
k = kr/r;
misc.omega = k*c_f;
misc.BC = 'SHBC';
msh.degree = 2;

pml.t = (R-r)/4;         % thickness of PML
pml.dirichlet = true;	% use homogeneous Dirichlet condition at Gamma_b (as opposed to homogeneous Neumann condition)
misc.r_a = R; 
pml.refinement = @(M) 2^(M-1)-1;
varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 3*2^(M-1)-1];
% varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1];
ffp.alpha_s = pi; % This is incorrectly set to pi in the paper
ffp.beta_s = 0;
ffp.r = R;
% ffp.alpha = linspace(0,2*pi,1000);
ffp.alpha = linspace(0,2*pi,361);
ffp.beta = 0;

para.plotResultsInParaview = 0;
ffp.calculateFarFieldPattern = 1;
ffp.plotFarField = false;     % If false, plots the near field instead
err.calculateVolumeError = 0;
err.calculateSurfaceError = 0;
loopParameters = {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','misc.omega'};

pml.sigmaType = [2,5,3];   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
pml.n = [2,2,1];
pml.sigmaType = 5;  % sigmaType = 5: sigma(xi) = gamma*xi^n
pml.n = 2;
pml.alpha = 30;
collectIntoTasks

misc.coreMethod = {'hp_FEM'};
varCol{1}.refinement = @(M) floor([(2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree/4-1]);
% collectIntoTasks

%% Plot results in paraview
misc.coreMethod = {'IGA'};
msh.explodeNURBS = true;   % Create patches from all C^0 interfaces
varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1, 2^(M-1)-1];
kr = 3*pi;
k = kr/r;
misc.omega = k*c_f;

pml.sigmaType = 3;  
pml.n = 1;
ffp.alpha_s = pi;
postPlot(1).plotResults  	= 0;
postPlot(1).printResults  	= 0;
postPlot(2).plotResults  	= 0;
postPlot(2).printResults  	= 0;
ffp.calculateFarFieldPattern = 0;
para.plotResultsInParaview = 1;
para.plotSubsets             = {'xy'};
para.plotFullDomain          = false;

para.extraXiPts              = 'round(2^(5-M)-1)';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(2^(5-M)-1)';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(2^(5-M)-1)';  % Extra visualization points in the zeta-direction per element
% collectIntoTasks



%% Plot frequency sweep
msh.explodeNURBS = false;   % Create patches from all C^0 interfaces
misc.scatteringCase = 'Sweep';
loopParameters = {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','ffp.alpha'};
noFreqs = 100;
% noFreqs = 2;
f_max = 2e3;
f = linspace(f_max/noFreqs,f_max,noFreqs);
misc.omega = 2*pi*f;
ffp.alpha_s = pi;
ffp.alpha = [0,pi];
% ffp.alpha = 0;
ffp.r = R;
pml.sigmaType = [2,5,3];   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
pml.n = [2,2,1];
% pml.sigmaType = 3;  
% pml.n = 1;
para.plotResultsInParaview = 0;
ffp.calculateFarFieldPattern = 1;
postPlot(3) = postPlot(1);
postPlot(3).plotResults  	= 1;
postPlot(3).printResults 	= 0;
postPlot(3).xname       	= 'misc.f';
postPlot(3).xScale          = 1;
postPlot(3).legendEntries 	= {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','ffp.alpha'};

postPlot(4)	= postPlot(3);
postPlot(4).yname        	= 'error_p';
postPlot(4).axisType        = 'semilogy';
% collectIntoTasks

misc.coreMethod = {'hp_FEM'};
varCol{1}.refinement = @(M) floor([(2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree/4-1]);
% collectIntoTasks



function addCommands_()

Mi_data = importdata('miscellaneous/refSolutions/Mi2021ilc.csv');
plot(Mi_data(:,1),Mi_data(:,2),'DisplayName','Mi2021ilc, M=2, $kr = \pi$');
hold on
plot(Mi_data(:,1),Mi_data(:,3),'DisplayName','Mi2021ilc, M=2, $kr = 2\pi$');
plot(Mi_data(:,1),Mi_data(:,4),'DisplayName','Mi2021ilc, M=2, $kr = 3\pi$');
plot(Mi_data(:,1),Mi_data(:,5),'DisplayName','Mi2021ilc, M=3, $kr = \pi$');
plot(Mi_data(:,1),Mi_data(:,6),'DisplayName','Mi2021ilc, M=3, $kr = 2\pi$');
plot(Mi_data(:,1),Mi_data(:,7),'DisplayName','Mi2021ilc, M=3, $kr = 3\pi$');
plot(Mi_data(:,1),Mi_data(:,8),'DisplayName','Mi2021ilc, M=4, $kr = \pi$');
plot(Mi_data(:,1),Mi_data(:,9),'DisplayName','Mi2021ilc, M=4, $kr = 2\pi$');
plot(Mi_data(:,1),Mi_data(:,10),'DisplayName','Mi2021ilc, M=4, $kr = 3\pi$');

legend('off');
legend('show','Interpreter','latex');

function addCommands2_()

analyticData = readLaTeXFormat('miscellaneous/refSolutions/Mi2021ilc_analytic.txt');
Mi_data = importdata('miscellaneous/refSolutions/Mi2021ilc.csv');

errors = 100*abs(Mi_data(:,2:end)-analyticData(:,2))./max(analyticData(:,2));
semilogy(Mi_data(:,1),errors(:,1),'DisplayName','Mi2021ilc, M=2, $kr = \pi$');
semilogy(Mi_data(:,1),errors(:,4),'DisplayName','Mi2021ilc, M=3, $kr = \pi$');
semilogy(Mi_data(:,1),errors(:,7),'DisplayName','Mi2021ilc, M=4, $kr = \pi$');
semilogy(Mi_data(:,1),errors(:,2),'DisplayName','Mi2021ilc, M=2, $kr = 2\pi$');
semilogy(Mi_data(:,1),errors(:,5),'DisplayName','Mi2021ilc, M=3, $kr = 2\pi$');
semilogy(Mi_data(:,1),errors(:,8),'DisplayName','Mi2021ilc, M=4, $kr = 2\pi$');
semilogy(Mi_data(:,1),errors(:,3),'DisplayName','Mi2021ilc, M=2, $kr = 3\pi$');
semilogy(Mi_data(:,1),errors(:,6),'DisplayName','Mi2021ilc, M=3, $kr = 3\pi$');
semilogy(Mi_data(:,1),errors(:,9),'DisplayName','Mi2021ilc, M=4, $kr = 3\pi$');

legend('off');
legend('show','Interpreter','latex');

