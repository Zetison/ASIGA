function m = knotRepetitions(Xi,xi)

m = zeros(size(xi));
for i = 1:numel(xi)
    m(i) = numel(find(Xi == xi(i)));
end