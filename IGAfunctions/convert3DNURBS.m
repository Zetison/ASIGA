function [Xi, Eta, Zeta, ...
          p, q, r, ...
          n, m, l, ...
          weights, controlPts] = convert3DNURBS(solid)

% convert GeoPDEs 3D NURBS object to the data structures
% used in IGA code.
% VP Nguyen, 2012.

p          = solid.degree(1);
q          = solid.degree(2);
r          = solid.degree(3);

Xi      = solid.knots{1};
Eta      = solid.knots{2};
Zeta      = solid.knots{3};

n     = length(Xi)-p-1;
m     = length(Eta)-q-1;
l     = length(Zeta)-r-1;

weights    = reshape(solid.coeffs(4,:,:),n*m*l,1);

controlPts = zeros(n*m*l,3);

count = 0;
for k=1:l
    for j=1:m
        controlPts(n*count+1:n*(count+1),:) = solid.coeffs(1:3,:,j,k)';
        count = count+1;
    end
end