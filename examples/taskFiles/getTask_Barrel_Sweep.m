function studies = getTask_Barrel_Sweep()

counter = 1;
studies = cell(0,1);
getDefaultTaskValues
noCoresToUse = 16;
saveStudies = false;

misc.model = 'Barrel_Sweep';

varCol = setBarrelParameters(1);
R = varCol{1}.R;
L = varCol{1}.L;
msh.meshFile = 'createNURBSmesh_Barrel';
msh.Xi = [0,0,0,1,1,2,2,3,3,3]/3;
pml.refinement = @(M) round((2^(M-1)-1)*pml.t/(R*2*pi/3));
ffp.calculateFarFieldPattern = false;
ffp.alpha_s = pi;
ffp.beta_s = 0;

ffp.beta = 0;
ffp.alpha = pi;

misc.solveForPtot = false;
misc.checkNURBSweightsCompatibility = false; 
warning('off','NURBS:weights')

prePlot.plot3Dgeometry = 1;
% prePlot.resolution = [20,20,20];
prePlot.resolution = [100,100,0];
prePlot.elementBasedSamples = 0;
prePlot.axis = 'off';
prePlot.view = [42,10];
prePlot.plotParmDir = 0;
prePlot.plotNormalVectors = 0;
prePlot.plotControlPolygon = 0;
prePlot.abortAfterPlotting = 1;                % Abort simulation after pre plotting
prePlot.coarseLinearSampling = prePlot.plotParmDir;
prePlot.plotFullDomain       = 1;

runTasksInParallel = ~prePlot.plot3Dgeometry;

msh.parm = 2;
msh.refineThetaOnly = false;

postPlot(1).xname       	= 'varCol{1}.kL';
postPlot(1).yname        	= 'surfaceError';
postPlot(1).plotResults  	= true;
postPlot(1).printResults 	= true;
postPlot(1).axisType        = 'semilogy';
postPlot(1).lineStyle   	= '-';
postPlot(1).fileDataHeaderX	= [];
postPlot(1).noXLoopPrms   	= 1;
postPlot(1).xLoopName     	= 'misc.omega';
postPlot(1).xScale          = 1;

msh.degree = 2;
misc.applyLoad = 'pointPulsation';
misc.method = {'BEM'};
misc.formulation = {'CCBIE','CHBIE','CBM','CCBIEC'};
% misc.formulation = {'CCBIE','CHBIE'};
misc.applyLoad = 'pointPulsation';
err.calculateSurfaceError = strcmp(misc.applyLoad,'pointPulsation');
misc.solveForPtot = ~strcmp(misc.applyLoad,'pointPulsation');
misc.BC = 'NBC';
msh.M = 4;
misc.scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering
% noPts = 100;
noPts = 1000;
kL_max = 15;
kL = linspace(kL_max/noPts,kL_max,noPts);
% Zeros of besselfunctions are found by utils/findBesselZeros.m
load('miscellaneous/besselZeros/besselJZeros.mat')
n3 = 1:100;
besselJZeros = reshape(besselJZeros,1,size(besselJZeros,1),size(besselJZeros,2));
eigenValuesCBIE = sqrt((besselJZeros/R).^2 + (n3.'*pi/L).^2);
eigenValuesCBIE = sort(unique(eigenValuesCBIE(:)));
eigenValuesCBIE(eigenValuesCBIE(:) > kL_max/L) = [];
load('miscellaneous/besselZeros/dbesselJZeros.mat')
n3 = 0:100;
noFuncs = size(dbesselJZeros,1);
dbesselJZeros = [0, dbesselJZeros(1,1:end-1); dbesselJZeros(2,:); [zeros(noFuncs-2,1), dbesselJZeros(3:end,1:end-1)]];
dbesselJZeros = reshape(dbesselJZeros,1,noFuncs,size(dbesselJZeros,2));
eigenValuesHBIE = sqrt((dbesselJZeros/R).^2 + (n3.'*pi/L).^2);
eigenValuesHBIE = sort(unique(eigenValuesHBIE(:)));
eigenValuesHBIE(eigenValuesHBIE(:) > kL_max/L) = [];
eigenValues = sort(unique([eigenValuesCBIE; eigenValuesHBIE]));
k = kL/L;
delta = 10/noPts*3;
% delta = 1e-2*kL_max/L;
k = [k, eigenValues.'];
for i = 1:numel(eigenValues)
    k = [k, eigenValues(i)+linspace(-delta/2,delta/2,round(noPts/10))];
%     k = [k, eigenValues(i)+linspace(-delta/2,delta/2,2)];
end
k(k <= 0) = [];
k(k > kL_max/L) = [];
k = sort(unique(k));

misc.omega = k*varCol{1}.c_f;
loopParameters = {'msh.M','msh.parm','misc.method','misc.formulation','misc.omega'};
collectIntoTasks


misc.method = {'BA'};
misc.formulation = {'SL2E'};
collectIntoTasks
