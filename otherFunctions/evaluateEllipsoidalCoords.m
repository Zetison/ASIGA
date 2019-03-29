function [r, theta, phi, T, c_1, c_2] = evaluateEllipsoidalCoords(x,y,z,Upsilon)



if abs(z/r) <= 1
    theta = acos(z/r);
elseif z/r > 1
    theta = acos(1);
else 
    theta = acos(-1);
end
phi = atan2(y,x);

if phi < 0
    phi = phi + 2*pi;
end