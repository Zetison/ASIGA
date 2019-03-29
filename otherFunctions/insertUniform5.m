function Y = insertUniform5(X, n)
% Y = insertUniform(X, n) inserts n uniform points distributed between X(1)
% and X(end) and stores it in Y (which does not include the values in X)

X = unique(X);
if n == 0
    Y = [];
    return
end
Y = zeros(n, 1);

n_arr = zeros(length(X)-1,1);
err_arr = zeros(length(X)-1,1);
for i = 1:length(n_arr)
    x = X(i+1)-X(i);
    n_arr(i) = floor(x*n);
    err_arr(i) = x/(n_arr(i)+1);
end
[~, I] = sort(err_arr,'descend');
for i = 1:n-sum(n_arr)
    n_arr(I(i)) = n_arr(I(i)) + 1;
end
idx = 0;
for i = 1:length(X)-1
    Y(idx+1:idx+n_arr(i)) = linspace2(X(i),X(i+1),n_arr(i))';
    idx = idx + n_arr(i);
end