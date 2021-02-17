function varCol = collectVariables(varCol,task)

varCol{1}.saveName       = task.saveName;
varCol{1}.scatteringCase = task.scatteringCase;
varCol{1}.M              = task.M;
varCol{1}.formulation    = task.formulation;
varCol{1}.coreMethod     = task.coreMethod;
varCol{1}.applyLoad      = task.applyLoad;
if isfield(task,'parm')
    varCol{1}.parm = task.parm;
end
if isfield(task,'internalPts')
    varCol{1}.internalPts = task.internalPts;
end
method = task.method;
if strcmp(method,'IENSG') || strcmp(method,'IE')
    varCol{1}.N = task.N;
    varCol{1}.IEbasis = task.IEbasis;
end
if strcmp(method,'ABC')
    varCol{1}.N = task.N;
end
  
for i = 1:numel(varCol)
    varCol{i}.extraGP = task.extraGP;
    varCol{i}.extraGPBEM = task.extraGPBEM;
    varCol{i}.buildMassMatrix = true;
    varCol{i}.buildStiffnessMatrix = ~strcmp(task.method,'BA');
    varCol{i}.applyBodyLoading = false;
    varCol{i}.extraGP = task.extraGP;
    varCol{i}.coreMethod = task.coreMethod;
    varCol{i}.formulation = task.formulation;
    switch varCol{i}.media
        case 'fluid'
            varCol{i}.dimension = 1;   
            varCol{i}.agpBEM = task.agpBEM;
            varCol{i}.quadMethodBEM = task.quadMethodBEM;
            varCol{i}.exteriorProblem = false;
            varCol{i}.model = task.model;
            varCol{i}.BC = task.BC;
            varCol{i}.useNeumanProj = task.useNeumanProj;
            varCol{i}.colBEM_C0 = task.colBEM_C0;
            varCol{i}.colMethod = task.colMethod;
            
            varCol{i}.solveForPtot = task.solveForPtot;
            if i == 1 && strcmp(task.method,'PML')
                varCol{i}.operator = 'PML';
            end
            varCol{i}.operator = 'Laplace';
            varCol{i}.fieldDimension = 1;
        case 'solid'
            varCol{i}.dimension = 3;   
            
            varCol{i}.solveForPtot = false;
            varCol{i}.operator = 'linearElasticity';
            varCol{i}.fieldDimension = 3;
    end  
    varCol{i} = convertNURBS(varCol{i});
end
varCol{1}.exteriorProblem = task.exteriorProblem;