scatteringCase = task.scatteringCase;
model = task.model;
method = task.method;
M = task.M;
f = task.f;

applyLoad = task.applyLoad;

alpha_s  = task.alpha_s;
beta_s  = task.beta_s;
alpha = task.alpha;
beta = task.beta;
r = task.r;
N_max = task.N_max;

degree = task.degree;

BC = task.BC;
plot2Dgeometry          = task.plot2Dgeometry;
plot3Dgeometry          = task.plot3Dgeometry;
calculateSurfaceError   = task.calculateSurfaceError;
calculateVolumeError    = task.calculateVolumeError;
calculateFarFieldPattern = task.calculateFarFieldPattern;
storeSolution           = task.storeSolution;
plotResultsInParaview   = task.plotResultsInParaview;
plotTimeOscillation     = task.plotTimeOscillation;
plotMesh                = task.plotMesh;
plotFarField            = task.plotFarField;
formulation             = task.formulation;
coreMethod              = task.coreMethod;
clearGlobalMatrices     = task.clearGlobalMatrices;
computeCondNumber       = task.computeCondNumber;
storeFullVarCol         = task.storeFullVarCol;
extraGP                 = task.extraGP;
extraGPBEM              = task.extraGPBEM;
agpBEM                  = task.agpBEM;
useROM                  = task.useROM;
exteriorProblem         = task.exteriorProblem;
varCol.initMeshFactXi  	= task.initMeshFactXi;
varCol.initMeshFactZeta	= task.initMeshFactZeta;
varCol.useNeumanProj    = task.useNeumanProj;
varCol.colBEM_C0        = task.colBEM_C0;
varCol.colMethod        = task.colMethod;
varCol.solveForPtot     = task.solveForPtot;
plotResidualError       = task.plotResidualError;
if isfield(task,'delta')
    varCol.delta = task.delta;
end

if isfield(task,'noVecsArr')
    noVecsArr = task.noVecsArr;
    noVecs = max(noVecsArr);
end
if isfield(task,'k_ROM')
    k_ROM = task.k_ROM;
end
if strcmp(method,'IENSG') || strcmp(method,'IE')
    if ~isfield(task,'N')
        switch formulation
            case {'PGC','BGC'}
                task.N = degree;
            case {'PGU','BGU'}
                task.N = 3;
        end
    end
    if ~isfield(task,'IEbasis')
        task.IEbasis = 'Chebyshev';
    end
else
    if ~isfield(task,'N')
        task.N = 2;
    end
end

if isfield(task,'parm')
    parm = task.parm;
end
if isfield(task,'internalPts')
    internalPts = task.internalPts;
end
    

if strcmp(method, 'BEM') && ~(strcmp(formulation, 'CCBIE') || strcmp(formulation, 'GCBIE') || ...
                              strcmp(formulation, 'CHBIE') || strcmp(formulation, 'GHBIE') ||  ...
                              strcmp(formulation, 'CBM') || strcmp(formulation, 'GBM') || ...
                              strcmp(formulation, 'CRCBIE1') || strcmp(formulation, 'GRCBIE1') || ...
                              strcmp(formulation, 'CRCBIE2') || strcmp(formulation, 'GRCBIE2') || ...
                              strcmp(formulation, 'CRCBIE3') || strcmp(formulation, 'GRCBIE3') || ...
                              strcmp(formulation, 'CCBIEC') || strcmp(formulation, 'GCBIEC') || ...
                              strcmp(formulation, 'CHBIEC') || strcmp(formulation, 'GHBIEC') ||  ...
                              strcmp(formulation, 'CBMC') || strcmp(formulation, 'GBMC') || ...
                              strcmp(formulation, 'CRCBIE1C') || strcmp(formulation, 'GRCBIE1C') || ...
                              strcmp(formulation, 'CRCBIE2C') || strcmp(formulation, 'GRCBIE2C') || ...
                              strcmp(formulation, 'CRCBIE3C') || strcmp(formulation, 'GRCBIE3C'))
	error(['The only "formulation"s available for "method = BEM" are CCBIE, CHBIE, CBM, GCBIE, GHBIE, GBM, '...
                'CRCBIE1, CRCBIE2, CRCBIE3, GRCBIE1, GRCBIE2 and GRCBIE3'])
end
if strcmp(method, 'IE') && ~(strcmp(formulation, 'PGU') || ...
                             strcmp(formulation, 'BGU') || ...
                             strcmp(formulation, 'PGC') || ...
                             strcmp(formulation, 'BGC'))
	error('The only "formulation"s available for "method = IE" are PGU, BGU, PGC and BGC')
end
if strcmp(method, 'IENSG') && ~(strcmp(formulation, 'PGU') || ...
                                strcmp(formulation, 'BGU') || ...
                                strcmp(formulation, 'PGC') || ...
                                strcmp(formulation, 'BGC'))
	error('The only "formulation"s available for "method = IENSG" are PGU, BGU, PGC and BGC')
end


varCol.method = method;





