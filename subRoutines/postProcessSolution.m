function [task,U] = postProcessSolution(task,UU)

if isfield(task.varCol{1},'noCols_tot')
    noCols_tot = task.varCol{1}.noCols_tot;
else
    noCols_tot = task.varCol{1}.noDofs;
end
if isfield(task.varCol{1},'allDofsToRemove')
    allDofsToRemove = task.varCol{1}.allDofsToRemove;
else
    allDofsToRemove = task.varCol{1}.dofsToRemove;
end
U = zeros(noCols_tot,size(UU,2));
U(setdiff(1:noCols_tot, allDofsToRemove'),:) = UU;   

h_max = -Inf;
for i = 1:numel(task.varCol)
    if isfield(task.varCol{1},'Aindices')
        task.varCol{i}.U = U(task.varCol{1}.Aindices{i,2},:);
    else
        task.varCol{i}.U = U;
    end
    h_max_i = findMaxElementDiameter(task.varCol{i}.patches);
    if h_max < h_max_i
        h_max = h_max_i;
    end
end
if ~isempty(allDofsToRemove)
    task.varCol = addSolutionToRemovedNodes(task.varCol);
end

if isfield(task.varCol{1},'c_f')
    task.h_max = h_max;
    task.nepw = task.varCol{1}.lambda./h_max;
    task.surfDofs = getNoSurfDofs(task);
    if task.varCol{1}.boundaryMethod
        task.tau = computeTau(task);
    end
end
for field = {'a','R','L'}
    if isfield(task.varCol{1},field{1})
        task.varCol{1}.(['k' field{1}]) = task.varCol{1}.(field{1})*task.varCol{1}.k;
    end
end