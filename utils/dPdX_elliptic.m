function A = dPdX_elliptic(f,r,theta)

a = r^2-f^2;
b = r^2-f^2*cos(theta)^2;

A = [a*cos(theta)/b      r*sqrt(a)*sin(theta)/b;
     -r*sin(theta)/b     sqrt(a)*cos(theta)/b];