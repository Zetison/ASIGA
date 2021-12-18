function studies = getTask_unitTest()


counter = 1;
studies = cell(0,1);
getDefaultTaskValues

misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
misc.model = 'S1';
varCol = setS1Parameters();
varCol = varCol(1);
noCoresToUse = 4;
msh.explodeNURBS = false;   % Create patches from all C^0 interfaces

msh.meshFile = 'createNURBSmesh_EL';
msh.parm = 1;
misc.checkNURBSweightsCompatibility = 0;
prePlot.plotGeometryInfo    = 0;       % Plot domain boundaries (i.e. Gamma, Gamma_a, Neumann, Dirichlet, ...)
err.calculateVolumeError    = true;

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

postPlot = [];

msh.meshFile = 'createNURBSmesh_EL';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
% msh.Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
msh.refineThetaOnly = 1;

connectedParameters = {};
% misc.method = {'IENSG','IE'};
misc.method = {'IE'};
misc.formulation = {'BGU'};
misc.coreMethod = {'IGA'};

msh.M = 5; % 4
misc.omega = 1000;
misc.BC = 'SHBC';
msh.degree = 2;
iem.boundaryMethod = 0;
iem.N = 64;

% misc.r_a = varCol{1}.R_i;
misc.r_a = varCol{1}.R_i*(1+2*pi/(32-pi));

varCol{1}.refinement = @(M) [2^(M-1)-1, 2^(M-1)-1, 2^(M-41)-1, 2^(M-4)-1];
% varCol{1}.refinement = @(M) [2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1, 2^(M+1)-1];
ffp.alpha_s = 0;
ffp.beta_s = pi/2;

ffp.calculateFarFieldPattern = 0;
err.calculateVolumeError = 1;
loopParameters = {'msh.M','misc.method'};

collectIntoTasks



