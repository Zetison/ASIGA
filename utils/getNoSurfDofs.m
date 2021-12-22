function surfDofs = getNoSurfDofs(task)

if task.varCol{1}.boundaryMethod
    switch task.misc.method
        case 'IENSG'
            surfDofs = task.dofs/task.iem.N;
        otherwise
            surfDofs = task.dofs;
    end
else
    varColBdry = meshBoundary(task.varCol{1},'Gamma');
    surfDofs = varColBdry.noSurfDofs - numel(varColBdry.dofsToRemove);
end