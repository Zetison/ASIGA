function studies = getTask_Mi2021ilc()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 19 to 22 in Mi2021ilc
% Mi2021ilc is available at https://doi.org/10.1016/j.cma.2021.113925

counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'Mi2021ilc';
noCoresToUse = 4;
msh.explodeNURBS = true;   % Create patches from all C^0 interfaces

msh.meshFile = 'createNURBSmesh_EL';
msh.parm = 1;
misc.checkNURBSweightsCompatibility = 0;
prePlot.plotGeometryInfo    = 0;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
err.calculateVolumeError    = true;
% misc.method = 'BA';
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% for BC = {'SHBC', 'SSBC','NNBC'}

prePlot.plotFullDomain   = false;        % Plot volumetric domains
prePlot.view             = [0,90];
prePlot.plotSubsets      = {'xy'};
% prePlot.plotSubsets      = {};
prePlot.plot3Dgeometry = 0;
prePlot.plot2Dgeometry = 0;
prePlot.plotControlPolygon = 0;       % Plot the control polygon for the NURBS mesh
prePlot.abortAfterPlotting = true;       % Abort simulation after pre plotting
% prePlot.colorFun = @(v) abs(norm2(v)-(r_a+t_PML));
prePlot.resolution = [20,20,0];
warning('off','NURBS:weights')

% postPlot(1).xname       	= 'nepw';
postPlot(1).xname       	= 'ffp.alpha';
postPlot(1).yname        	= 'abs_p';
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 0;
postPlot(1).axisType        = 'plot';
postPlot(1).lineStyle   	= '-';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 0;
postPlot(1).xScale       	= 180/pi;
postPlot(1).legendEntries 	= {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','varCol{1}.k'};

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
msh.M = 4; % 4
varCol{1}.media = 'fluid'; % Media; % solid or fluid (Navier equation or Helmholtz equation)
varCol{1}.R_i = r;
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
varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1, 2^(M-1)-1];
% varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1];
ffp.alpha_s = 0; % This is incorrectly set to pi in the paper
ffp.beta_s = 0;
ffp.r = R;
ffp.alpha = linspace(0,2*pi,1000);
ffp.beta = 0;

para.plotResultsInParaview = 0;
ffp.calculateFarFieldPattern = 1;
ffp.plotFarField = false;     % If false, plots the near field instead
err.calculateVolumeError = 0;
err.calculateSurfaceError = 0;
loopParameters = {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','misc.omega'};

pml.sigmaType = [2,3];   	% sigmaType = 1: sigma(xi) = xi*exp(gamma*xi), sigmaType = 2: sigma(xi) = C*xi^n, sigmaType = 3: sigma(xi) = C/(1-xi)^n
pml.n = [2,1];
% pml.sigmaType = 3;  
% pml.n = 1;
collectIntoTasks

misc.coreMethod = {'hp_FEM'};
varCol{1}.refinement = @(M) floor([(2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree/4-1]);
collectIntoTasks

%% Plot results in paraview
misc.coreMethod = {'IGA'};
varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1, 2^(M-1)-1];
kr = 3*pi;
k = kr/r;
misc.omega = k*c_f;

ffp.alpha_s = pi;
postPlot(1).plotResults  	= 0;
postPlot(1).printResults  	= 0;
ffp.calculateFarFieldPattern = 0;
para.plotResultsInParaview = 1;
para.plotSubsets             = {'xy'};
para.plotFullDomain          = false;

para.extraXiPts              = 'round(2^(5-M)-1)';  % Extra visualization points in the xi-direction per element
para.extraEtaPts             = 'round(2^(5-M)-1)';  % Extra visualization points in the eta-direction per element
para.extraZetaPts            = 'round(2^(5-M)-1)';  % Extra visualization points in the zeta-direction per element
collectIntoTasks



%% Plot frequency sweep
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
para.plotResultsInParaview = 0;
ffp.calculateFarFieldPattern = 1;
postPlot(1).plotResults  	= 1;
postPlot(1).printResults 	= 0;
postPlot(1).xname       	= 'misc.f';
postPlot(1).xScale          = 1;
postPlot(1).legendEntries 	= {'msh.M','msh.degree','pml.sigmaType','misc.method','misc.coreMethod','ffp.alpha'};
collectIntoTasks

misc.coreMethod = {'hp_FEM'};
varCol{1}.refinement = @(M) floor([(2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree-1, (2^(M+1)-1+msh.degree)/msh.degree/4-1]);
collectIntoTasks



