function I = simpson(f,a,b,m)

h = (b-a)/m;
sum1 = 0;
for j = 1:m/2-1
    x = a + (2*j)*h;
    sum1 = sum1 + f(x);
end
sum2 = 0;
for j = 1:m/2
    x = a + (2*j-1)*h;
    sum2 = sum2 + f(x);
end

I = h/3*(f(a) + 2*sum1 + 4*sum2 + f(b));