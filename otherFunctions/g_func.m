function g = g_func(a,x,y)

g = zeros(size(x));
[I, J] = size(a);
for i = 1:I
    m = i-1;
    for j = 1:J
        n = j-1;
        g = g + a(i,j)*x.^m.*y.^n;
    end
end