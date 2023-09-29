function Y = insertUniform2(X, n)
% Y = insertUniform(X, n) inserts n uniform points between each unique
% value in X and stores it in Y (which does not include the values in X)

X = unique(X);
if n < 1
    Y = [];
    return
end
Y = zeros(n*(length(X)-1), 1);

for i = 1:length(X)-1
    idx = n*(i-1);
    Y(idx+1:idx+n) = linspace2(X(i),X(i+1),n)';
end