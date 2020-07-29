
varCol{1}.saveName = saveName;
varCol{1}.scatteringCase = scatteringCase;
varCol{1}.M = M;
varCol{1}.formulation = formulation;
varCol{1}.coreMethod = coreMethod;
varCol{1}.applyLoad = applyLoad;
varCol{1}.alpha_s = alpha_s;
varCol{1}.beta_s = beta_s;
varCol{1}.alpha = alpha;
varCol{1}.beta = beta;
varCol{1}.analyticSolutionExist = analyticSolutionExist;
if isfield(task,'parm')
    varCol{1}.parm = task.parm;
end
if isfield(task,'internalPts')
    varCol{1}.internalPts = task.internalPts;
end

if strcmp(method,'IENSG') || strcmp(method,'IE')
    varCol{1}.N = task.N;
    varCol{1}.IEbasis = task.IEbasis;
    varCol{1} = generateCoeffMatrix(varCol{1});
end
if strcmp(method,'ABC')
    varCol{1}.N = task.N;
end
  
for i = 1:numel(varCol)
    varCol{i}.extraGP = extraGP;
    varCol{i}.extraGPBEM = extraGPBEM;
    switch varCol{i}.media
        case 'fluid'
            varCol{i}.dimension = 1;   
            varCol{i}.agpBEM = agpBEM;
            varCol{i}.quadMethodBEM =task.quadMethodBEM;
            varCol{i}.extraGP = extraGP;
            varCol{i}.coreMethod = coreMethod;
            varCol{i}.formulation = formulation;
            varCol{i}.exteriorProblem = false;
            varCol{i}.model = model;
            varCol{i}.BC = BC;
            varCol{i}.useNeumanProj = useNeumanProj;
            varCol{i}.colBEM_C0 = colBEM_C0;
            varCol{i}.colMethod = colMethod;
            varCol{i}.solveForPtot = solveForPtot;
            varCol{i}.useSolidDomain = useSolidDomain;
        case 'solid'
            varCol{i}.dimension = 3;   
            varCol{i}.extraGP = extraGP;
            varCol{i}.solveForPtot = false;
            varCol{i}.coreMethod = coreMethod;
            varCol{i}.formulation = formulation;
    end  
    varCol{i} = convertNURBS(varCol{i});
end
varCol{1}.exteriorProblem = exteriorProblem;