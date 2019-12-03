
task.model = model;
task.plot2Dgeometry = plot2Dgeometry;
task.plot3Dgeometry = plot3Dgeometry;
task.calculateSurfaceError = calculateSurfaceError;
if calculateSurfaceError
    if exist('LpOrder','var')
        task.LpOrder = LpOrder;
    else
        task.LpOrder = 2;
    end
end
if exist('internalPts','var')
    task.internalPts = internalPts;
end
task.calculateVolumeError = calculateVolumeError;
task.calculateSurfEnrgErr = calculateSurfEnrgErr;
task.calculateFarFieldPattern = calculateFarFieldPattern;
task.useROM = useROM;
if useROM
    storeSolution = true;
end
task.storeSolution = storeSolution;
task.plotResultsInParaview = plotResultsInParaview;
task.plotTimeOscillation = plotTimeOscillation;
task.plotMesh = plotMesh;
task.plotFarField = plotFarField;
task.M = M;
task.r = r;
task.f = f;
task.applyLoad = applyLoad;
task.model = model;
task.scatteringCase = scatteringCase;
task.alpha_s = alpha_s;
task.beta_s = beta_s;
task.formulation = formulation;
task.storeFullVarCol = storeFullVarCol;
task.clearGlobalMatrices = clearGlobalMatrices;
task.N_max = N_max;
task.extraGP = extraGP;
task.extraGPBEM = extraGPBEM;
task.agpBEM = agpBEM;
task.exteriorProblem = exteriorProblem;
task.initMeshFactXi = initMeshFactXi;
task.initMeshFactZeta = initMeshFactZeta;
task.method = method;
task.useNeumanProj = useNeumanProj;
task.colBEM_C0 = colBEM_C0;
task.plotResidualError = plotResidualError;
task.colMethod = colMethod;
task.solveForPtot = solveForPtot;
if exist('r_a','var')
    task.r_a = r_a;
end
if exist('N','var')
    task.N = N;
end
if exist('noVecsArr','var')
    task.noVecsArr = noVecsArr;
end
if exist('k_ROM','var')
    task.k_ROM = k_ROM;
end
if exist('useDGP','var')
    task.useDGP = useDGP;
else
    task.useDGP = false;
end
if exist('IEbasis','var')
    task.IEbasis = IEbasis;
end

task.coreMethod = coreMethod;
task.degree = degree;
task.BC = BC;
task.alpha = alpha;
task.beta = beta;
if exist('parm','var')
    task.parm = parm;
end
if exist('delta','var')
    task.delta = delta;
end
task.quadMethodBEM = quadMethodBEM;

task.computeCondNumber = computeCondNumber;
if ~exist('loopParameters','var')
    loopParameters = {};
end
if (~isnan(alpha_s(1)) || ~isnan(beta_s(1))) && (strcmp(scatteringCase,'MS') || strcmp(scatteringCase,'Sweep'))
    error(['For monostatic scattering alpha_s and beta_s should not be given (they should be defined through alpha and beta). ' ...
          'Note that alpha_s = alpha and beta_s = beta.'])
end
loopParametersArr = cell(length(loopParameters),1);
for i = 1:length(loopParameters)
    loopParametersArr{i} = task.(loopParameters{i});
    if ischar(loopParametersArr{i})
        error('Strings in loopParametersArr must be stored in cell arrays.');
    end
end

studies(counter).loopParameters = loopParameters;
studies(counter).loopParametersArr = loopParametersArr;
studies(counter).runTasksInParallel = runTasksInParallel;
if exist('basisROMcell','var')
    studies(counter).basisROMcell = basisROMcell;
    studies(counter).k_ROM = k_ROM;
    studies(counter).noVecsArr = noVecsArr;
end

studies(counter).tasks = createTasks([], 1, task, 1, loopParameters, loopParametersArr);

counter = counter + 1;

