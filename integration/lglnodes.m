function [x, w, P] = lglnodes(p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods. 
%
% Reference on LGL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N = p+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x = cos(pi*(p:-1:0)/p);

% The Legendre Vandermonde Matrix
P = zeros(N,N);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold = 2;

while max(abs(x-xold)) > eps
    xold = x;
        
    P(1,:) = 1;    
    P(2,:) = x;
    
    for k = 2:p
        P(k+1,:) = ((2*k-1)*x.*P(k,:) - (k-1)*P(k-1,:))/k;
    end
     
    x = xold - (x.*P(N,:) - P(p,:))./(N*P(N,:));
             
end

w = 2./(p*N*P(N,:).^2);