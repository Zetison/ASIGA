function r = getNACA2(s,t,s2)
% returns the n'th derivative of NACA profile with respect to xi = sqrt(x)

f = @(xi) getNACA(xi,t);
xi = invertNACA(s,s2,t);
r = [xi.^2, f(xi)];
