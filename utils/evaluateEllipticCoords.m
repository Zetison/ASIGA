function [r, theta] = evaluateEllipticCoords(x,y,f)

T_f = x^2+y^2+f^2;

r = 1/sqrt(2)*sqrt(T_f+sqrt(T_f^2 - 4*f^2*x^2));

theta = atan2(r*y,x*sqrt(r^2-f^2));
% theta = acos(x/r);
% if y < 0
%     theta = 2*pi-theta;
% end