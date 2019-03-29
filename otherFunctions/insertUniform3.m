function Y = insertUniform3(X, n, p)
% Y = insertUniform(X, n) inserts n uniform points between each unique
% value in X which is repeated p times and stores it in Y

uniqeX = unique(X);
distr = histc(X, uniqeX);

X_interp = uniqeX(distr >= p);
 
Y = zeros((n+1)*length(X_interp)-n, 1);

for i = 1:length(X_interp)-1
    idx = (n+1)*i-n;
    Y(idx) =  X_interp(i);
    Y(idx+1:idx+n) = linspace2(X_interp(i),X_interp(i+1),n)';
end
Y(end) = X_interp(end);