function [r, theta, phi, d1, d2] = evaluateProlateCoords(X,Upsilon)
x = X(:,1);
y = X(:,2);
z = X(:,3);
d1 = sqrt(x.^2+y.^2+(z+Upsilon).^2);
d2 = sqrt(x.^2+y.^2+(z-Upsilon).^2);

r = 0.5*(d1+d2);
theta = zeros(size(x));
indices = abs(z./r) <= 1;
theta(indices) = acos(z(indices)./r(indices));
indices = z./r > 1;
theta(indices) = 0;
indices = z./r < -1;
theta(indices) = pi;

phi = atan2(y,x);
indices = phi < 0;
phi(indices) = phi(indices) + 2*pi;  % require 0 <= phi <= 2*pi