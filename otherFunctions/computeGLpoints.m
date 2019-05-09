function [Q,W] = computeGLpoints(n)

% mp.Digits(100);

Q = cell(n,1);
W = cell(n,1);
for N = 1:n
    N
    [Q{N}, W{N}] = lgwt(N);
end

function [x, w] = lgwt(N)
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
d = mp.Digits;
N = mp(N);
N=N-1;
N1=N+1; N2=N+2;
xu=linspace(mp('-1'),mp('1'),N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2,'mp');
% Derivative of LGVM
Lp=zeros(N1,N2,'mp');
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=mp('2');
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0)) > 10^(-d)
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end
% Linear map from[-1,1] to [a,b]
x=(-(1-y)+(1+y))/2;      
% Compute the weights
w=2./((1-y.^2).*Lp.^2)*(N2/N1)^2;
