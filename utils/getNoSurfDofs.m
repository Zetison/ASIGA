function surfDofs = getNoSurfDofs(task)

if task.varCol{1}.boundaryMethod
    switch task.misc.method
        case 'IENSG'
            surfDofs = task.dofs/task.iem.N;
        otherwise
            surfDofs = task.dofs;
    end
else
    counter = 0;
    surfDofs = [];
    for patch = 1:task.varCol{1}.noPatches
        number = task.varCol{1}.patches{patch}.nurbs.number;
        surfDofsPatch = 1:prod(number(1:2));
        surfDofs = [surfDofs, surfDofsPatch+counter];
        counter = counter + prod(number);
    end
    surfDofs = numel(setdiff(surfDofs,task.varCol{1}.dofsToRemove));
end