function [x,w]=lgwt(N)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [-1,1] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [-1,1]
% which you can evaluate at any x in [-1,1]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; 
N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
x=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

x0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(x-x0)./abs(x0)) > eps
    L(:,1) = 1;
    
    L(:,2) = x;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*x.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-x.*L(:,N2) )./(1-x.^2);   
    
    x0=x;
    x=x0-L(:,N2)./Lp;
end

% Compute the weights
w=2./((1-x.^2).*Lp.^2)*(N2/N1)^2;



