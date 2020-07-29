function idxMap = findMatchIndices(X, epsilon)

idxMap = cell(1,0);
counter = 1;
foundMatch = false;
for i = 1:size(X,1)-1
    % Check whether or not a nod is already considered
    if ~isempty(find(cell2mat(idxMap) == i))
        continue;
    end
    for j = i+1:size(X,1)
        if norm(X(i,:) - X(j,:)) < epsilon
            if ~foundMatch
                idxMap{counter} = [i, j];
                foundMatch = true;
            else
                idxMap{counter} = [idxMap{counter} j];
            end
        end
    end
    if foundMatch
        counter = counter + 1;
        foundMatch = false;
    end
end