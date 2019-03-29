function A = vandermonde(v)

N = length(v);
A = zeros(N);

A(:,1) = ones(N,1);

for i = 2:N
    A(:,i) = v.^(i-1);
end