function varCol2 = extractVarColFields(task,varCol)


varCol2{1}.alpha = varCol{1}.alpha;
varCol2{1}.beta = varCol{1}.beta;
varCol2{1}.k = varCol{1}.k;
varCol2{1}.f = varCol{1}.f;
varCol2{1}.c_f = varCol{1}.c_f;
if ~strcmp(task.method,'RT') && ~strcmp(task.method,'KDT')
    varCol2{1}.surfDofs = varCol{1}.surfDofs;
    varCol2{1}.dofs = varCol{1}.dofs;
    varCol2{1}.dofsAlg = (varCol{1}.dofs)^(1/3);
    varCol2{1}.h_max = varCol{1}.h_max;
    if isfield(varCol{1},'tau')
        varCol2{1}.tau = varCol{1}.tau;
    end
    varCol2{1}.nepw = varCol{1}.nepw;
    varCol2{1}.tot_time = varCol{1}.tot_time;
    varCol2{1}.noElems = varCol{1}.noElems;
    varCol2{1}.totNoElems = varCol{1}.totNoElems;
    if isfield(varCol{1},'N')
        varCol2{1}.N = varCol{1}.N;
    end
    if isfield(varCol{1},'timeBuildSystem')
        varCol2{1}.timeBuildSystem = varCol{1}.timeBuildSystem;
    end
    if isfield(varCol{1},'timeSolveSystem')
        varCol2{1}.timeSolveSystem = varCol{1}.timeSolveSystem;
    end
    if isfield(varCol{1},'totNoQP')
        varCol2{1}.totNoQP = varCol{1}.totNoQP;
    end
    if isfield(varCol{1},'totNoQPnonPolar')
        varCol2{1}.totNoQPnonPolar = varCol{1}.totNoQPnonPolar;
    end
end 
varCol2{1}.analyticSolutionExist = varCol{1}.analyticSolutionExist;
varCol2{1}.boundaryMethod = varCol{1}.boundaryMethod;