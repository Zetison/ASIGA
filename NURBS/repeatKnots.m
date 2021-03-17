function task = repeatKnots(task)
for i_varCol = 1:numel(task.varCol) % assume coreMethod to be the same in all domains
    nurbsPatches = task.varCol{i_varCol}.nurbs;
    if ~iscell(nurbsPatches)
        nurbsPatches = {nurbsPatches};
    end
    for patch = 1:numel(nurbsPatches)
        nurbs = nurbsPatches(patch);
        switch task.misc.coreMethod
            case {'linear_FEM', 'h_FEM', 'hp_FEM', 'C0_IGA', 'SEM'}
                d_p = nurbs{1}.d_p;
                knots = nurbs{1}.knots;
                newKnots = cell(1,d_p);
                for j = 1:d_p  
                    for i = 1:length(knots{j})
                        mm = length(find(knots{j} == knots{j}(i)));  
                        newKnots{j} = [newKnots{j}, knots{j}(i)*ones(nurbs{1}.degree(j)-mm,1)];
                    end
                end
                nurbs = insertKnotsInNURBS(nurbs,newKnots);
        end
        nurbsPatches(patch) = nurbs;
    end
    task.varCol{i_varCol}.nurbs = nurbsPatches;
end