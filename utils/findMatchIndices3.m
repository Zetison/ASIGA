function idxMap = findMatchIndices3(X, epsilon,M,N,candidatesIndices)

idxMatrix = zeros(M,N);


ii = 1;
foundMatch = false;
for i = 1:length(candidatesIndices)
    cnd_i = candidatesIndices(i);
    % Check whether or not a node is already considered
    if ~isempty(find(idxMatrix == cnd_i))
        continue;
    end
    for j = i+1:length(candidatesIndices)
        cnd_j = candidatesIndices(j); 
        if (abs(X(cnd_i,1) - X(cnd_j,1)) < epsilon) && ...
           (abs(X(cnd_i,2) - X(cnd_j,2)) < epsilon) && ...
           (abs(X(cnd_i,3) - X(cnd_j,3)) < epsilon)
            if ~foundMatch
                idxMatrix(ii,1:2) = [cnd_i, cnd_j];
                foundMatch = true;
                jj = 3;
            else
                idxMatrix(ii,jj) = cnd_j;
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