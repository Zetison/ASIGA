function task = controlTask(task)

formulation = task.misc.formulation;
scatteringCase = task.misc.scatteringCase;
method = task.misc.method;
applyLoad = task.misc.applyLoad;

if task.ffp.calculateFarFieldPattern && task.msh.refineThetaOnly && max(task.ffp.extraGP) < 7
    warning('Consider using more integration points (ffp.extrGP) in the parametric direction not refined to ensure the fundamental functions (which are not necessarily axisymmetric) are integrated properly in the far field calculation routine')
end
if strcmp(method,'PML') && isnan(task.pml.t) 
    error('The PML thickness pml.t must be set')
end
if strcmp(method,'PML') && isnan(task.misc.r_a) 
    error('The distance to the PML layer must be set')
end
if strcmp(method,'PML') && ~isfield(task.pml,'refinement') 
    error('The refinement field must be specified using PML')
end
if strcmp(task.misc.method,'BEM') && ~task.misc.solveForPtot && ~strcmp(task.misc.BC,'NBC')...
        && ~strcmp(applyLoad,'radialPulsation')
    warning('It is reccomended to use solveForPtot = true for BEM')
end
if task.misc.solveForPtot && ~strcmp(applyLoad,'planeWave')
    error('p_inc does not solve the interior problem in the case applyLoad=radialPulsation. For this reason one must have solveForPtot=false here.')
end
if strcmp(scatteringCase,'MS') && ~strcmp(applyLoad,'planeWave')
    error('MS only applies for planeWave')
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
if (strcmp(method, 'IE') || strcmp(method, 'IENSG')) && ~(strcmp(formulation, 'PGU') || ...
                             strcmp(formulation, 'BGU') || ...
                             strcmp(formulation, 'PGC') || ...
                             strcmp(formulation, 'BGC') || ...
                             strcmp(formulation, 'WBGC') || ...
                             strcmp(formulation, 'WBGU'))
	error('The only "formulation"s available for "method = IE" or method = IENSG are PGU, BGU, PGC, BGC, WBGC and WBGU')
end
if strcmp(method, 'PML') && ~(strcmp(formulation, 'GSB') || strcmp(formulation, 'STD'))
	error('The only "formulation"s available for "method = PML" is GSB (STD is only implemented for testing purposes)')
end
if strcmp(method,'BA') && strcmp(scatteringCase,'MS') && (numel(task.ffp.alpha_s) > 1 || numel(task.ffp.beta_s) > 1)
    error('This is case is not implemented. The best approximation method must be combined with "scatteringCase = BI"')
end
if strcmp(method,'PML') && isnan(task.pml.gamma) && task.rom.useROM
    error('gamma is not set: gamma cannot be frequency dependent when using ROM')
end
task.misc.storeFullVarCol = task.misc.storeFullVarCol || task.rom.useROM;

if strcmp(task.misc.scatteringCase,'MS') && (~isnan(task.ffp.alpha_s(1)) || ~isnan(task.ffp.beta_s(1)))
    warning(['For monostatic scattering alpha_s and beta_s should not be given (they should be defined through alpha and beta). ' ...
           'The given variables will be overwritten as alpha_s = alpha and beta_s = beta.'])
end
if (isnan(task.ffp.alpha_s(1)) || isnan(task.ffp.beta_s(1))) && ~strcmp(task.misc.scatteringCase,'MS') && strcmp(applyLoad,'planeWave')
    error('Incident direction is not set: alpha_s = NaN and/or beta_s = NaN')
end
if task.misc.solveForPtot && ~(strcmp(task.misc.method,'BEM') || strcmp(task.misc.method,'BA'))
    error('solveForPtot can only be used with method = BEM or method = BA')
end
if (task.pml.sigmaType == 3 || task.pml.sigmaType == 4) && ~task.pml.dirichlet
    error('sigmaType == 3 and sigmaType == 4 requires Dirichlet conditions at Gamma_b')
end

if task.err.calculateVolumeError && isBoundaryMethod(task)
    warning('There is no volumetric errors to compute for boundary methods')
end


