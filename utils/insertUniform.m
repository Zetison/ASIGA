function Y = insertUniform(X, n)
% Y = insertUniform(X, n) inserts n uniform points between each unique
% value in X and stores it in Y

X = unique(X);
Y = zeros((n+1)*length(X)-n, 1);

for i = 1:length(X)-1
    idx = (n+1)*i-n;
    Y(idx) =  X(i);
    Y(idx+1:idx+n) = linspace2(X(i),X(i+1),n)';
end
Y(end) = X(end);