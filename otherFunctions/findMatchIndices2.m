function idxMap = findMatchIndices2(X, epsilon,M,N)

idxMatrix = zeros(M,N);


ii = 1;
foundMatch = false;
for i = 1:size(X,1)-1
    % Check whether or not a node is already considered
    if ~isempty(find(idxMatrix == i))
        continue;
    end
    for j = i+1:size(X,1)
        if (abs(X(i,1) - X(j,1)) < epsilon) && ...
           (abs(X(i,2) - X(j,2)) < epsilon) && ...
           (abs(X(i,3) - X(j,3)) < epsilon)
            if ~foundMatch
                idxMatrix(ii,1:2) = [i, j];
                foundMatch = true;
                jj = 3;
            else
                idxMatrix(ii,jj) = j;
                jj = jj + 1;
            end
        end
    end
    if foundMatch
        ii = ii + 1;
        foundMatch = false;
    end
end

noConnections = length(find(idxMatrix(:,1) ~= 0));

idxMatrix = idxMatrix(1:noConnections,:);

idxMap = cell(1,noConnections);
if noConnections == 0
    return
end

for ii = 1:noConnections   
    nonZeroElementIndices = (idxMatrix(ii,:) ~= 0);
    idxMap{ii} = idxMatrix(ii,nonZeroElementIndices);
end