function U = addSolutionToRemovedNodes_new(U, varCol, d)
if nargin < 3
    d = varCol.dimension;
end
gluedNodes = varCol.gluedNodes;

for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        for ii = 1:d
            U(d*(gluedNodes{i}(j)-1)+ii,:) = U(d*(parentIdx-1)+ii,:);
        end
    end
end