function varCol = addSolutionToRemovedNodes(varCol)

for patch = 1:numel(varCol)
    if isfield(varCol{patch},'U')
        gluedNodes = varCol{patch}.gluedNodes;
        d = varCol{patch}.dimension;
        for i = 1:length(gluedNodes)
            parentIdx = gluedNodes{i}(1);
            for j = 2:length(gluedNodes{i})
                for ii = 1:d
                    varCol{patch}.U(d*(gluedNodes{i}(j)-1)+ii,:) = varCol{patch}.U(d*(parentIdx-1)+ii,:);
                end
            end
        end
    end
end