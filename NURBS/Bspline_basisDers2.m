function [N, dNdxi, d2Ndxi2, d3Ndxi3] = Bspline_basisDers2(i, xi, p, Xi)
% This routine compute the p+1 nonzero basis functions and corresponding
% derivatives at xi

% Input
%       i:      knot span index corresonding to xi
%       p:      the degree of the B-Spline/NURBS
%       xi:     the value for which we want to evaluate the Bspline
%       Xi:     an open knot vector of size n+p+1

% Output
%       N:      array of the p+1 B-spline functions evaluated at xi  



computeDers = nargout > 1;
noxi = numel(xi);
N = zeros(noxi,p+1,class(xi));
N(:,1) = 1;
saved = ones(noxi,1,class(xi));

for j = 2:p+1
    % For k = 1 there is no dependence on N(k-1) of the previous run.
    for k = 1:j
        % Compute N_{i-j+k,j-1} according to the Cox-deBoor formula
        temp = zeros(noxi,1);
        if k ~= j
            temp = (Xi(i+k)-xi)/(Xi(i+k)-Xi(i-j+k+1)).*N(:,k);
        end
        if k ~= 1
            temp = temp + (xi-Xi(i-j+k))/(Xi(i+k-1)-Xi(i-j+k)).*saved;
        end
        saved = N(:,k);
        N(:,k) = temp;
    end
    if j == p
        N_tilde = N(:,1:end-1);
    elseif j == p-1
        N_tilde2 = N(:,1:end-2);
    elseif j == p-2
        N_tilde3 = N(:,1:end-3);
    end
end
if p == 3
    N_tilde3 = 1;
end

if p == 2
    N_tilde2 = 1;
end

if p == 1
    N_tilde = 1;
end

if computeDers
    %% Calculate first derivatives
    dNdxi = zeros(noxi,p+1,class(xi));
    a = 1./(Xi(i+1:i+p)-Xi(i-p+1:i));
    a = repmat(a,noxi,1);
    dNdxi(:,1:p)   = -N_tilde.*a;
    dNdxi(:,2:p+1) = dNdxi(:,2:p+1) + N_tilde.*a;
    dNdxi = p*dNdxi;

    %% Calculate second derivatives
    d2Ndxi2 = zeros(noxi,p+1,class(xi));
    if p > 1
        a2 = 1./(Xi(i+1:i+p-1)-Xi(i-p+2:i));
        a2 = repmat(a2,noxi,1);
        d2Ndxi2(:,1:p-1) =                    N_tilde2.*a(:,1:end-1).*a2;
        d2Ndxi2(:,2:p)   = d2Ndxi2(:,2:p)   - N_tilde2.*a2.*(a(:,1:end-1) + a(:,2:end));
        d2Ndxi2(:,3:p+1) = d2Ndxi2(:,3:p+1) + N_tilde2.*a(:,2:end).*a2;
        d2Ndxi2 = p*(p-1)*d2Ndxi2;
    end

    %% Calculate third derivatives
    if nargout == 4
        d3Ndxi3 = zeros(noxi,p+1,class(xi));
        if p > 2
            a3 = 1./(Xi(i+1:i+p-2)-Xi(i-p+3:i));
            a3 = repmat(a3,noxi,1);

            d3Ndxi3(:,1:p-2) =                  - N_tilde3.*a(:,1:end-2).*a2(:,1:end-1).*a3; % 4
            d3Ndxi3(:,2:p-1) = d3Ndxi3(:,2:p-1) + N_tilde3.*a3.*(a(:,2:end-1).*a2(:,2:end) + (a(:,1:end-2) + a(:,2:end-1)).*a2(:,1:end-1)); % 3
            d3Ndxi3(:,3:p)   = d3Ndxi3(:,3:p)   - N_tilde3.*a3.*((a(:,2:end-1) + a(:,3:end)).*a2(:,2:end) + a(:,2:end-1).*a2(:,1:end-1)); % 2
            d3Ndxi3(:,4:p+1) = d3Ndxi3(:,4:p+1) + N_tilde3.*a(:,3:end).*a2(:,2:end).*a3; % 1
            d3Ndxi3 = p*(p-1)*(p-2)*d3Ndxi3;
        end
    end
end





