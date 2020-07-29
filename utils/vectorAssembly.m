function F = vectorAssembly(values,idx,N)
[n, m, l] = size(values);
F = zeros(N,l);
for j = 1:m
    for i = 1:n
        temp = values(i,j,:);
        F(idx(i,j),:) = F(idx(i,j),:) + temp(:).';
    end
end