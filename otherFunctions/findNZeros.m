function x_roots = findNZeros(f,N,dx)

x = 0.1;
while isnan(f(x)) || isinf(f(x))
    x = x+dx;
end
x_roots = zeros(N,1);
counter = 1;
noZerosFound = 0;
f_prev = f(x);
while noZerosFound < N
    x = x + dx;
    f_val = f(x);
    if sign(f_prev) ~= sign(f_val)
        x_roots(counter) = bisection(f, x-dx, x, floor(log(dx/10^(-13))/log(2)), x*10^(-13));
%         x_roots(counter) = fzero(f,x - dx/2);
        x = x_roots(counter) + dx;
        counter = counter + 1;
        noZerosFound = noZerosFound + 1;
    end
    f_prev = f_val;
end