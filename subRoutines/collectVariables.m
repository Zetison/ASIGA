function task = collectVariables(task)
  
for i = 1:numel(task.varCol)
    task.varCol{i}.buildMassMatrix = true;
    task.varCol{i}.buildStiffnessMatrix = ~strcmp(task.misc.method,'BA');
    task.varCol{i}.applyBodyLoading = false;
    switch task.varCol{i}.media
        case 'fluid'
            task.varCol{i}.dimension = 1;   
            task.varCol{i}.operator = 'Laplace';
            task.varCol{i}.fieldDimension = 1;
        case 'solid'
            task.varCol{i}.dimension = 3;   
            task.varCol{i}.solveForPtot = false;
            task.varCol{i}.operator = 'linearElasticity';
            task.varCol{i}.fieldDimension = 3;
    end  
    task.varCol{i} = convertNURBS(task.varCol{i});
end