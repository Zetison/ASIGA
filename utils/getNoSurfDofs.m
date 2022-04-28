function surfDofs = getNoSurfDofs(task)

if task.varCol{1}.boundaryMethod
    surfDofs = task.FEdofs;
else
    varColBdry = meshBoundary(task.varCol{1},'Gamma');
    surfDofs = varColBdry.noSurfDofs - numel(varColBdry.dofsToRemove);
end