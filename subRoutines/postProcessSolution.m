function task = postProcessSolution(task)

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
noCols_tot_red = size(task.UU,1);
task.UU = [task.UU; zeros(noCols_tot-noCols_tot_red,size(task.UU,2))];
task.UU(setdiff(1:noCols_tot, allDofsToRemove'),:) = task.UU(1:noCols_tot_red,:);
task.UU(allDofsToRemove,:) = 0;

for i = 1:numel(task.varCol)
    if isfield(task.varCol{1},'Aindices')
        task.varCol{i}.U = task.UU(task.varCol{1}.Aindices{i,2},:);
    else
        task.varCol{i}.U = task.UU;
    end
end
task = rmfield(task,'UU');
if ~isempty(allDofsToRemove)
    task.varCol = addSolutionToRemovedNodes(task.varCol);
end