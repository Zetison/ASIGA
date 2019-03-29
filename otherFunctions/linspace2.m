function x = linspace2(a,b,n)
dx = (b-a)/(n+1);
x = zeros(1,n);
for i = 1:n
    x(i) = a + i*dx;
end