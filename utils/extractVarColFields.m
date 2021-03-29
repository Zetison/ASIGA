function varCol2 = extractVarColFields(task,varCol)


varCol2{1}.alpha = task.ffp.alpha;
varCol2{1}.beta = task.ffp.beta;
if isfield(task.ffp,'theta')
    varCol2{1}.theta = task.ffp.theta;
end
varCol2{1}.f = task.misc.omega/(2*pi);
varCol2{1}.k = task.misc.omega/varCol{1}.c_f;
if isfield(task,'dofs')
    varCol2{1}.dofsAlg = (task.dofs)^(1/3);
end

varCol2{1}.totNoElems = task.totNoElems;    
fieldCell = {'noElems', 'tot_time', 'nepw', 'tau', 'h_max', 'surfDofs', 'dofs', 'N', 'timeBuildSystem', 'timeSolveSystem','totNoQP','totNoQPnonPolar','boundaryMethod'};

for field = fieldCell
    if isfield(varCol{1},field{1})
        varCol2{1}.(field{1}) = varCol{1}.(field{1});
    end
end