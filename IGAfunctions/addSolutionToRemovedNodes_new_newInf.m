function U = addSolutionToRemovedNodes_new_newInf(U, varCol)

gluedNodes = varCol.gluedNodes;
noDofs = varCol.noDofs;
N = varCol.N;

for m = 1:N
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            U(gluedNodes{i}(j)+(m-1)*noDofs) = U(parentIdx+(m-1)*noDofs);
        end
    end
end