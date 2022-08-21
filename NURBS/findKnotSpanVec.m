function mid = findKnotSpanVec(n, p, xi, Xi)
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
npts = numel(xi);
Xi = Xi.'; % Make Xi a column vector

low = p*ones(npts,1);
high = (n + 1)*ones(npts,1);
endIndices = xi == Xi(end);
mid = floor((low + high) / 2);
mid(endIndices) = n;
if all(endIndices)
    return
end
indices = and(or(xi < Xi(mid), xi >= Xi(mid+1)),~endIndices);
indices_high = false(npts,1);
indices_low = false(npts,1);
while any(indices)
    indices_high(~indices) = false;
    indices_high(indices) = xi(indices) < Xi(mid(indices));
    high(indices_high) = mid(indices_high);

    indices_low(~indices) = false;
    indices_low(indices) = ~indices_high(indices);
    low(indices_low) = mid(indices_low);

    mid(indices) = floor((low(indices) + high(indices)) / 2);
    indices = and(or(xi < Xi(mid), xi >= Xi(mid+1)),~endIndices);
end