function [varCol,U] = postProcessSolution(varCol,UU)

noCols_tot = varCol{1}.noCols_tot;
U = zeros(noCols_tot,size(UU,2));
U(setdiff(1:noCols_tot, varCol{1}.allDofsToRemove'),:) = UU;   
h_max = -Inf;
for i = 1:numel(varCol)
    varCol{i}.U = U(varCol{1}.Aindices{i,2},:);
    h_max_i = findMaxElementDiameter(varCol{i}.patches);
    if h_max < h_max_i
        h_max = h_max_i;
    end
end
varCol = addSolutionToRemovedNodes(varCol);

varCol{1}.h_max = h_max;
varCol{1}.nepw = varCol{1}.lambda./h_max;
varCol{1}.surfDofs = getNoSurfDofs(varCol{1});
if varCol{1}.boundaryMethod
    varCol{1}.tau = computeTau(varCol{1});
end