function task = computeDerivedQuantities(task)


totNoElems = 0;
FEdofs = 0;
h_max = -Inf;
for i = 1:task.noDomains
    totNoElems = totNoElems + task.varCol{i}.noElems;
    FEdofs = FEdofs + task.varCol{i}.noDofs - numel(task.varCol{i}.dofsToRemove);
    h_max_i = findMaxElementDiameter(task.varCol{i}.patches);
    if h_max < h_max_i
        h_max = h_max_i;
    end
end
task.h_max = h_max;
task.totNoElems = totNoElems;
task.FEdofs = FEdofs;

if isfield(task.varCol{1},'c_f')
    task.nepw = task.varCol{1}.lambda./task.h_max;
    task.surfDofs = getNoSurfDofs(task);
    if task.varCol{1}.boundaryMethod
        task.tau = computeTau(task);
    end
end
for field = {'a','R','L'}
    if isfield(task.varCol{1},field{1})
        task.(['k' field{1}]) = task.varCol{1}.(field{1})*task.varCol{1}.k;
    end
end
if isfield(task.misc,'omega')
    task.f = task.misc.omega/(2*pi);
    task.k = task.misc.omega/task.varCol{1}.c_f;
end
if isfield(task,'dofs')
    task.dofsAlg = (task.dofs)^(1/3);
end