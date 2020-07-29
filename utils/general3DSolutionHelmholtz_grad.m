function grad_u = general3DSolutionHelmholtz_grad(x,y,z,k, A, B, C, D, E, F)
% Solution with notation from 
% http://mathworld.wolfram.com/HelmholtzDifferentialEquationCartesianCoordinates.html
grad_u = zeros(3,1);

[L, M] = size(E);

for l = 1:L
    for m = 1:M
        lambda = sqrt(k^2+l^2+m^2);
        grad_u(1) = grad_u(1) + (A(l)*l*exp(l*x) - B(l)*l*exp(-l*x)).*(C(m)*exp(m*y) ...
            + D(m)*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
        grad_u(2) = grad_u(2) + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*m*exp(m*y) ...
            - D(m)*m*exp(-m*y)).*(E(l,m)*exp(-1i*lambda*z) + F(l,m)*exp(1i*lambda*z));
        grad_u(3) = grad_u(3) + (A(l)*exp(l*x) + B(l)*exp(-l*x)).*(C(m)*exp(m*y) ...
            + D(m)*exp(-m*y)).*(E(l,m)*-1i*lambda*exp(-1i*lambda*z) + F(l,m)*1i*lambda*exp(1i*lambda*z));
    end
end