function task = extractVarColFields(task)


varCol{1}.alpha = task.ffp.alpha;
varCol{1}.beta = task.ffp.beta;
if isfield(task.ffp,'theta')
    varCol{1}.theta = task.ffp.theta;
end
varCol{1}.f = task.misc.omega/(2*pi);
varCol{1}.k = task.misc.omega/task.varCol{1}.c_f;
if isfield(task,'dofs')
    varCol{1}.dofsAlg = (task.dofs)^(1/3);
end

varCol{1}.totNoElems = task.totNoElems;    
fieldCell = {'noElems', 'tot_time', 'nepw', 'tau', 'h_max', 'surfDofs', 'dofs', 'N', 'timeBuildSystem', 'timeSolveSystem','totNoQP','totNoQPnonPolar','boundaryMethod','kL','ka','kR'};

for field = fieldCell
    if isfield(task.varCol{1},field{1})
        varCol{1}.(field{1}) = task.varCol{1}.(field{1});
    end
end
task.varCol = varCol;