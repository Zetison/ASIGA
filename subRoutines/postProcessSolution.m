function [varCol,U] = postProcessSolution(varCol,UU)

if isfield(varCol{1},'noCols_tot')
    noCols_tot = varCol{1}.noCols_tot;
else
    noCols_tot = varCol{1}.noDofs;
end
if isfield(varCol{1},'allDofsToRemove')
    allDofsToRemove = varCol{1}.allDofsToRemove;
else
    allDofsToRemove = varCol{1}.dofsToRemove;
end
U = zeros(noCols_tot,size(UU,2));
U(setdiff(1:noCols_tot, allDofsToRemove'),:) = UU;   

h_max = -Inf;
for i = 1:numel(varCol)
    if isfield(varCol{1},'Aindices')
        varCol{i}.U = U(varCol{1}.Aindices{i,2},:);
    end
    h_max_i = findMaxElementDiameter(varCol{i}.patches);
    if h_max < h_max_i
        h_max = h_max_i;
    end
end
varCol = addSolutionToRemovedNodes(varCol);

if isfield(varCol{1},'lambda')
    varCol{1}.h_max = h_max;
    varCol{1}.nepw = varCol{1}.lambda./h_max;
    varCol{1}.surfDofs = getNoSurfDofs(varCol{1});
    if varCol{1}.boundaryMethod
        varCol{1}.tau = computeTau(varCol{1});
    end
end