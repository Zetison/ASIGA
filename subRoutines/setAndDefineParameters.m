
   
switch model
    case {'SS', 'PS', 'S1', 'S1_P2', 'S3', 'S5', 'S13', 'S15', 'S35', 'IL'}
        analyticSolutionExist = true;
    otherwise
        analyticSolutionExist = false;
end
if strcmp(model(end-1:end),'_P')
    analyticSolutionExist = true;
end
if strcmp(applyLoad,'radialPulsation') || strcmp(applyLoad,'SimpsonTorus')
    analyticSolutionExist = true;
end    
if (calculateSurfaceError || calculateVolumeError) && ~analyticSolutionExist
    error('The errors cannot be computed without an analytic solution')
end


if isfield(task,'parm')
    parms = setParameters(model,task.parm);
    varCol.parm = parm;
else
    parms = setParameters(model, []);
end
if isfield(task,'r_a')
    varCol.r_a = task.r_a;
end

P_inc = parms.P_inc;
rho_s = parms.rho_s;
rho_f = parms.rho_f;
c_f = parms.c_f;
E = parms.E;
nu = parms.nu;

isSphericalShell = any(strcmp(model, {'SS_P', 'SS', 'S1', 'S1_P', 'S1_P2', 'S3', 'S5', 'IL'}));
varCol.isSphericalShell = isSphericalShell;
if isSphericalShell
    parms.R_i = parms.R_o - parms.t; % Inner radius of shell
    R_i = parms.R_i; 
end
switch BC
    case {'SHBC','NBC'}
        E = E(1:end-1);
        rho_s = rho_s(1:end-1);
        rho_f = rho_f(1:end-1);
        c_f = c_f(1:end-1);
        nu = nu(1:end-1);
        if isSphericalShell
            R_i = R_i(1:end-1);
        end
        useSolidDomain = false;
        useInnerFluidDomain = false;
    case 'SFBC'
        rho_f = rho_f(1:end-1);
        c_f = c_f(1:end-1);
        if isSphericalShell
            R_i = R_i(1:end-1);
        end
        useSolidDomain = true;
        useInnerFluidDomain = false;
    case 'SSBC'
        rho_f = rho_f(1:end-1);
        c_f = c_f(1:end-1);
        useSolidDomain = true;
        useInnerFluidDomain = false;
    case 'NNBC'
        useSolidDomain = true;
        useInnerFluidDomain = true;
    otherwise
        useSolidDomain = false;
        useInnerFluidDomain = false;        
end
if useSolidDomain
    C = zeros(6,6);
    C(1:3,1:3) = E/(1+nu)/(1-2*nu)*...
                 [1-nu nu nu;
                  nu 1-nu nu; 
                  nu nu 1-nu];
    C(4:6,4:6) = 0.5*E/(1+nu)*eye(3);
    varCol_solid.C = C;
end
if useInnerFluidDomain
    if ~useSolidDomain
        error('An intermediate solid domain must be used when considering a seperate inner fluid domain!')
    end
end

switch method
    case {'IE','ABC'}
        boundaryMethod = false;
    case {'IENSG','BEM','KDT','MFS','RT'}
        boundaryMethod = true;
    case 'BA'
        switch formulation
            case {'SL2E','SL2Etot'}
                boundaryMethod = true;
            case {'VL2E','VL2Etot'}
                boundaryMethod = false;
            otherwise
                error('Not implemented')
        end
end

if boundaryMethod && useSolidDomain
    error('Not yet implemented!')
end
varCol.boundaryMethod = boundaryMethod;

switch scatteringCase
    case 'MS'
        alpha_s = alpha;
        beta_s = beta;
end

varCol.useSolidDomain = useSolidDomain;
varCol.useInnerFluidDomain = useInnerFluidDomain;
varCol.coreMethod = coreMethod;
