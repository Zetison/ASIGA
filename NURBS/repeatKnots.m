function task = repeatKnots(task)
switch task.misc.coreMethod
    case {'linear_FEM', 'h_FEM', 'hp_FEM', 'C0_IGA', 'SEM'}
        for i_varCol = 1:numel(task.varCol) % assume coreMethod to be the same in all domains
            task.varCol{i_varCol}.nurbs = repeatNURBSknots(task.varCol{i_varCol}.nurbs);
        end
end