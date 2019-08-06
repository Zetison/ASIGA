function surfDofs = getNoSurfDofs(varCol)

if varCol.boundaryMethod
    surfDofs = varCol.dofs;
else
    counter = 0;
    surfDofs = [];
    for patch = 1:varCol.noPatches
        number = varCol.patches{patch}.nurbs.number;
        surfDofsPatch = 1:prod(number(1:2));
        surfDofs = [surfDofs, surfDofsPatch+counter];
        counter = counter + prod(number);
    end
    surfDofs = numel(setdiff(surfDofs,varCol.dofsToRemove));
end