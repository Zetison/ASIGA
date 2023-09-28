function m = knotRepetitions(Xi,xi,Eps)

if nargin < 3
    Eps = 1e-10;
end
m = zeros(size(xi));
for i = 1:numel(xi)
    m(i) = numel(find(abs(Xi - xi(i)) < Eps));
end