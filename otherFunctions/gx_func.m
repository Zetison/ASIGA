function g = gx_func(a,x,y)

g = zeros(size(x));
[I, J] = size(a);
for i = 2:I
    m = i-1;
    for j = 1:J
        n = j-1;
        g = g + a(i,j)*m*x.^(m-1).*y.^n;
    end
end