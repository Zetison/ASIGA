function g = gy_func(a,x,y)

g = zeros(size(x));
[I, J] = size(a);
for i = 1:I
    m = i-1;
    for j = 2:J
        n = j-1;
        g = g + a(i,j)*n*x.^m.*y.^(n-1);
    end
end