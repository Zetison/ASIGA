function [N, dNdxi] = Bspline_basisDers(i, xi, p, Xi)
error('Depricated. Use BsplineBasis instead')
% This routine compute the p+1 nonzero basis functions and corresponding
% derivatives at xi

% Input
%       i:      knot span index corresonding to xi
%       p:      the degree of the B-Spline/NURBS
%       xi:     the value for which we want to evaluate the Bspline
%       Xi:     an open knot vector of size n+p+1

% Output
%       N:      array of the p+1 B-spline functions evaluated at xi  




N = zeros(1,p+1);
N(1) = 1;
saved = 1;

for j = 2:p+1
    % For k = 1 there is no dependence on N(k-1) of the previous run.
    for k = 1:j
        % Compute N_{i-j+k,j-1} according to the Cox-deBoor formula
        temp = 0;
        if k ~= j
            temp = (Xi(i+k)-xi)/(Xi(i+k)-Xi(i-j+k+1))*N(k);
        end
        if k ~= 1
            temp = temp + (xi-Xi(i-j+k))/(Xi(i+k-1)-Xi(i-j+k))*saved;
        end
        saved = N(k);
        N(k) = temp;
    end
    if j == p
        N_tilde = N(1:end-1);
    end
end
if p == 1
    N_tilde = 1;
end

dNdxi = zeros(1,p+1);
dNdxi(1:p)   = -p*N_tilde./(Xi(i+1:i+p)-Xi(i-p+1:i));
dNdxi(2:p+1) = dNdxi(2:p+1) + p*N_tilde./(Xi(i+1:i+p)-Xi(i+1-p:i));


