function g = dtu(x)

g = zeros(size(x));
indices = and(0 < x, x < 1);
x = x(indices);

f1 = exp(-1./x);
f2 = exp(-1./(1-x));
df1 = f1./x.^2;
df2 = f2./(1-x).^2;
g(indices) = (df1.*f2+f1.*df2)./(f1+f2).^2;

return

x = linspace(-1,2,1000);
plot(x,tu(x),x,dtu(x));
