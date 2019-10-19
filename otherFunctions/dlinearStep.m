function f = dlinearStep(x)

f = zeros(size(x));
f(and(-eps <= x, x <= 1+eps)) = 1;
return

x = linspace(-1,2,1000);
plot(x,linearStep(x),x,dlinearStep(x));