function u = general3DSolutionHelmholtz(x,y,z,k, A, B, C, D, E, F)
% Solution with notation from 
% http://mathworld.wolfram.com/HelmholtzDifferentialEquationCartesianCoordinates.html
u = 0;

[L, M] = size(E);

for l = 1:L
    for m = 1:M
        lambda = sqrt(k^2+l^2+m^2);
        u = u + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) ...
            + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
    end
end