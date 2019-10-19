function f = linearStep(x)

f = zeros(size(x));
f(x >= 1) = 1;
indices = and(0 < x, x < 1);
x = x(indices);
f(indices) = x;
return

x = linspace(-1,2,1000);
plot(x,linearStep(x));