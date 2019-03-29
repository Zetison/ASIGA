function A = dXdP(r,theta,phi,f)

c = sqrt(r^2-f^2);

A = [r*sin(theta)*cos(phi)/c         r*sin(theta)*sin(phi)/c    cos(theta);
     c*cos(theta)*cos(phi)           c*cos(theta)*sin(phi)      -r*sin(theta);
     -c*sin(theta)*sin(phi)          c*sin(theta)*cos(phi)      0];