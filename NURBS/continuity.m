function k = continuity(Xi,xi,p)

k = Inf(size(xi));
for i = 1:numel(xi)
    m = numel(find(Xi == xi(i)));
    if m > 0
        k(i) = p-m;
    end
end