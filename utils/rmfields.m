function varCol = rmfields(varCol,fieldCell)

for i = 1:numel(fieldCell)
    if iscell(varCol)
        for j = 1:numel(varCol)
            if isfield(varCol{j},fieldCell{i})
                varCol{j} = rmfield(varCol{j},fieldCell{i});
            end
        end
    else
        if isfield(varCol,fieldCell{i})
            varCol = rmfield(varCol,fieldCell{i});
        end
    end
end