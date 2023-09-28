function varCol = addSolutionToRemovedNodes(varCol)

for i_domain = 1:numel(varCol)
    if isfield(varCol{i_domain},'U')
        gluedNodes = varCol{i_domain}.gluedNodes;
        d = varCol{i_domain}.dimension;
        for i = 1:length(gluedNodes)
            parentIdx = gluedNodes{i}(1);
            for j = 2:length(gluedNodes{i})
                for ii = 1:d
                    varCol{i_domain}.U(d*(gluedNodes{i}(j)-1)+ii,:) = varCol{i_domain}.U(d*(parentIdx-1)+ii,:);
                end
            end
        end
    end
end