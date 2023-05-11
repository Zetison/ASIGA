function Y = insertNonUniform(X,N)
% Y = insertUniform(X, n) inserts N(i) uniform points between [X(i),X(i+1)]  and stores it in Y (which does not include the values in X)

X = unique(X);
if isempty(N)
    Y = [];
    return
end
if numel(N) ~= numel(X)-1
    error('The length of N must be equal to the intervals/elements of X')
end
Y = zeros(sum(N), 1);
idx = 0;
for i = 1:length(X)-1
    Y(idx+1:idx+N(i)) = linspace2(X(i),X(i+1),N(i))';
    idx = idx + N(i);
end