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

degree = task.degree;

BC                      = task.BC;
prePlot                 = task.prePlot;
para                    = task.para;
calculateSurfaceError   = task.calculateSurfaceError;
calculateVolumeError    = task.calculateVolumeError;
calculateFarFieldPattern = task.calculateFarFieldPattern;
storeSolution           = task.storeSolution;
plotFarField            = task.plotFarField;
formulation             = task.formulation;
coreMethod              = task.coreMethod;
clearGlobalMatrices     = task.clearGlobalMatrices;
computeCondNumber       = task.computeCondNumber;
storeFullvarCol         = task.storeFullVarCol;
extraGP                 = task.extraGP;
extraGPBEM              = task.extraGPBEM;
agpBEM                  = task.agpBEM;
useROM                  = task.useROM;
exteriorProblem         = task.exteriorProblem;
useNeumanProj           = task.useNeumanProj;
colBEM_C0               = task.colBEM_C0;
colMethod               = task.colMethod;
solveForPtot            = task.solveForPtot;
varCol                  = task.varCol;
varCol{1}.initMeshFactXi   = task.initMeshFactXi;
varCol{1}.initMeshFactZeta = task.initMeshFactZeta;
varCol{1}.refineThetaOnly  = task.refineThetaOnly;
varCol{1}.parm          = task.parm;
varCol{1}.r_a           = task.r_a;
varCol{1}.N             = task.N;
varCol{1}.method        = method;
varCol{1}.model         = model;
varCol{1}.useROM        = useROM;
varCol{1}.N_max         = task.N_max;
varCol{1}.ie_Zeta       = task.ie_Zeta;
varCol{1}.IElocSup      = task.IElocSup;
varCol{1}.p_ie          = task.p_ie;
for i = 1:numel(varCol)
    varCol{i}.progressBars = task.progressBars;
end

plotResidualError       = task.plotResidualError;
if isfield(task,'delta')
    varCol{1}.delta = task.delta;
end
if isfield(task,'Xi')
    varCol{1}.Xi = task.Xi;
end
if isfield(task,'c_x')
    varCol{1}.c_x = task.c_x;
end
if isfield(task,'c_y')
    varCol{1}.c_y = task.c_y;
end
if isfield(task,'c_z')
    varCol{1}.c_z = task.c_z;
end
if isfield(task,'P_inc')
    varCol{1}.P_inc = task.P_inc;
else
    varCol{1}.P_inc = 1;
end
if strcmp(task.method,'BEM') && ~task.solveForPtot && ~strcmp(task.BC,'NBC')...
        && ~strcmp(task.applyLoad,'radialPulsation')
    warning('It is reccomended to use solveForPtot = true for BEM')
end
if task.solveForPtot && ~strcmp(task.applyLoad,'planeWave')
    error('p_inc does not solve the interior problem in the case applyLoad=radialPulsation. For this reason one must have solveForPtot=false here.')
end
if isfield(task,'noVecsArr')
    noVecsArr = task.noVecsArr;
    noVecs = max(noVecsArr);
end
if isfield(task,'k_ROM')
    k_ROM = task.k_ROM;
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
if strcmp(method,'BA') && strcmp(scatteringCase,'MS')
    error('This is case is not implemented. The best approximation method must be combined with "scatteringCase = BI"')
end






