function i = findKnotSpan(n, p, xi, Xi)
% This routine finds the knot span corresponding to a given value xi.
% The method uses a sequencial search algorithm

% Input
%       n:    the number of control points
%       p:    the degree of the B-Spline
%       xi:   the value for which we want to find the knot span
%       Xi:   an open knot vector of size n+p+1

% Output
%       i: index of knot such that xi is an element in [Xi(i) Xi(i+1))

% Check for xi = Xi(end)
% if xi < 0 || xi > 1
%     keyboard
% end
xi = real(xi);
if xi == Xi(end)
    i = n;
    return;
end

low = p;
high = n + 1;
mid = floor((low + high) / 2);
while xi < Xi(mid) || xi >= Xi(mid+1)
    if xi < Xi(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high) / 2);
end

i = mid;