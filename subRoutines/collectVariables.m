
omega = 2*pi*f;
k = (1./c_f.')*omega;
lambda = 2*pi./k;  % Wave length for outer fluid domain

varCol.saveName = saveName;
varCol.scatteringCase = scatteringCase;
varCol.BC = BC;
varCol.model = model;
varCol.M = M;
varCol.f = f;
varCol.rho_f = rho_f(1);
varCol.c_f = c_f(1);
varCol.k = k(1,:);
varCol.formulation = formulation;
varCol.coreMethod = coreMethod;
varCol.omega = omega;
varCol.applyLoad = applyLoad;
varCol.alpha_s = alpha_s;
varCol.beta_s = beta_s;
varCol.alpha = alpha;
varCol.beta = beta;
varCol.extraGP = extraGP;
varCol.extraGPBEM = extraGPBEM;
varCol.agpBEM = agpBEM;
varCol.parms = parms;
varCol.P_inc = P_inc;
varCol.analyticSolutionExist = analyticSolutionExist;
varCol.exteriorProblem = exteriorProblem;
if isfield(task,'parm')
    varCol.parm = task.parm;
end
if isfield(task,'internalPts')
    varCol.internalPts = task.internalPts;
end

if strcmp(method,'IENSG') || strcmp(method,'IE')
    varCol.N = task.N;
    varCol.IEbasis = task.IEbasis;
    varCol = generateCoeffMatrix(varCol);
end
if strcmp(method,'ABC')
    varCol.N = task.N;
end

%% Convert NURBS data    
varCol.dimension = 1;    
varCol = convertNURBS(fluid, varCol);

if useSolidDomain 
    varCol_solid.dimension = 3;     
    varCol_solid = convertNURBS(solid, varCol_solid);   
    if varCol_solid.patches{1}.nurbs.degree(1) ~= varCol.patches{1}.nurbs.degree(1) || ...
       varCol_solid.patches{1}.nurbs.degree(2) ~= varCol.patches{1}.nurbs.degree(2)
        error('This case has not been implemented')
    end
    varCol_solid.omega = omega;
    varCol_solid.rho_s = rho_s;
end
if useInnerFluidDomain     
    varCol_fluid_i.dimension = 1;   
    varCol_fluid_i = convertNURBS(fluid_i, varCol_fluid_i);
    varCol_fluid_i.omega = omega;
    varCol_fluid_i.k = k(:,2);
    varCol_fluid_i.rho_f = rho_f(2);
end
