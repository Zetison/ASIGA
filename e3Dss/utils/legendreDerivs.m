function [P, dP, d2P] = legendreDerivs(n, x, P, dP, d2P)
% This function returns n'th Legendre polynomial evaluated at x It is
% assumed that the input arguments P, dP, d2P contains the evaluation of
% the Legendre polynomial, it's derivative and it's double derivative at x
% of order n-1 and n-2 stored at P(1,:) and P(2,:) respectively (and so
% on).

if n == 0
    P = [zeros(size(x),class(x)); ones(size(x),class(x))];
elseif n == 1
    P(1,:) = P(2,:);
    P(2,:) = x;
else
    temp = P(2,:);
    P(2,:) = ((2*n-1)*x.*P(2,:) - (n-1)*P(1,:))/n;
    P(1,:) = temp;
end
if nargin >= 4
    if n == 1
        dP = [zeros(size(x),class(x)); ones(size(x),class(x))];
    elseif n == 2
        dP(1,:) = dP(2,:);
        dP(2,:) = 3*x;
    elseif n > 2
        temp = dP(2,:);
        dP(2,:) = ((2*n-1)*(P(1,:) + x.*dP(2,:)) - (n-1)*dP(1,:))/n;
        dP(1,:) = temp;
    end
end
if nargin >= 5
    if n == 2
        d2P = [zeros(size(x),class(x)); 3*ones(size(x),class(x))];
    elseif n == 3
        d2P(1,:) = d2P(2,:);
        d2P(2,:) = 15*x;
    elseif n > 3
        temp = d2P(2,:);
        d2P(2,:) = ((2*n-1)*(2*dP(1,:) + x.*d2P(2,:)) - (n-1)*d2P(1,:))/n;
        d2P(1,:) = temp;
    end
end