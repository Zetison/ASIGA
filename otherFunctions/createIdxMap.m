function idxMap = createIdxMap(N)

idxMap = zeros(N,N-1);
for j = 1:N
    idxMap(j,:) = [1:j-1,j+1:N];
end