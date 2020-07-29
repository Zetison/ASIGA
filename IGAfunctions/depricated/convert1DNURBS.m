function [Xi, p, n, weights, controlPts] = convert1DNURBS(solid)
error('Depricated, use convertNURBS instead')

% convert GeoPDEs 3D NURBS object to the data structures
% used in IGA code.
% VP Nguyen, 2012.

p          = solid.degree(1);

Xi      = solid.knots;

n     = length(Xi)-p-1;

weights    = reshape(solid.coeffs(2,:),n,1);

controlPts = reshape(solid.coeffs(1,:),n,1);
