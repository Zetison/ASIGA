function A = dXdP(r,theta,phi,Upsilon)

c = sqrt(r^2-Upsilon^2);

A = [r*sin(theta)*cos(phi)/c	c*cos(theta)*cos(phi)	-c*sin(theta)*sin(phi);
     r*sin(theta)*sin(phi)/c	c*cos(theta)*sin(phi)  	c*sin(theta)*cos(phi);
     cos(theta)                 -r*sin(theta)           0];