function [r, theta, phi, d1, d2] = evaluateProlateCoords(x,y,z,Upsilon)

d1 = sqrt(x^2+y^2+(z+Upsilon)^2);
d2 = sqrt(x^2+y^2+(z-Upsilon)^2);

r = 0.5*(d1+d2);

if abs(z/r) <= 1
    theta = acos(z/r);
elseif z/r > 1
    theta = 0;
else 
    theta = pi;
end
phi = atan2(y,x);

if phi < 0 % require 0 <= phi <= 2*pi
    phi = phi + 2*pi;
end