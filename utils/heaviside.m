function h = heaviside(x)

h = zeros(size(x));
h(x >= -eps) = 1;
