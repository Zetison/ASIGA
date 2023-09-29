function ders = BsplineBasis(i, xi, p, Xi, n)
% This routine compute the p+1 nonzero basis functions and corresponding
% n derivatives at xi
% It is based on Algorithm A2.3 in The NURBS Book Piegl1997tnb available at
% https://doi.org/10.1016/S0167-8396(98)00028-4

% Input
%       i:      knot span index corresonding to xi
%       p:      the degree of the B-Spline
%       xi:     the value for which we want to evaluate the Bspline
%       Xi:     an open knot vector of size n+p+1
%       n:      number of derivatives to compute

% Output
%       ders: 	array of the p+1 B-spline functions and its derivatives 
%           	evaluated at xi  
if nargin < 5
    n = 0;
end
Xi = Xi.'; % Make Xi a column vector
noxi = numel(xi);
left = zeros(noxi,p+1);
right = zeros(noxi,p+1);
ndu = zeros(noxi,p+1,p+1);
a = zeros(noxi,p+1,2);
ndu(:,1,1) = ones(noxi,1);
for j = 1:p
    left(:,j+1) = xi-Xi(i-j+1);
    right(:,j+1) = Xi(i+j)-xi;
    saved = zeros(noxi,1);
    for r = 0:j-1
        ndu(:,r+1,j+1) = right(:,r+2)+left(:,j-r+1);
        temp = ndu(:,j,r+1)./ndu(:,r+1,j+1);
        ndu(:,j+1,r+1) = saved+right(:,r+2).*temp;
        saved = left(:,j-r+1).*temp;
    end
    ndu(:,j+1,j+1) = saved;
end
ders = zeros(noxi,p+1,n+1);
for j = 0:p
    ders(:,j+1,1) = ndu(:,p+1,j+1);
end
for r = 0:p
    s1 = 0; 
    s2 = 1;
    a(:,1,1) = ones(noxi,1);
    for k = 1:n
        d = zeros(noxi,1);
        rk = r-k;
        pk = p-k;
        if r >= k
            a(:,1,s2+1) = a(:,1,s1+1)./ndu(:,rk+1,pk+2);
            d = a(:,1,s2+1).*ndu(:,pk+1,rk+1);
        end
        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        if r-1 <= pk
            j2 = k-1;
        else
            j2 = p-r;
        end
        for j = j1:j2
            a(:,j+1,s2+1) = (a(:,j+1,s1+1)-a(:,j,s1+1))./ndu(:,rk+j+1,pk+2);
            d = d + a(:,j+1,s2+1).*ndu(:,pk+1,rk+j+1);
        end
        if r <= pk
            a(:,k+1,s2+1) = -a(:,k,s1+1)./ndu(:,r+1,pk+2);
            d = d + a(:,k+1,s2+1).*ndu(:,pk+1,r+1);
        end
        ders(:,r+1,k+1) = d;
        j = s1;
        s1 = s2;
        s2 = j;
    end
end
r = p;
for k = 1:n
    for j = 0:p
        ders(:,j+1,k+1) = ders(:,j+1,k+1)*r;
    end
    r = r*(p-k);
end