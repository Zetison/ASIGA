function A = dPdX(x,y,z,Upsilon,r,d1,d2)

A = [x*(d1+d2)/(2*d1*d2), y*(d1+d2)/(2*d1*d2), (z*(d1+d2)+Upsilon*(d2-d1))/(2*d1*d2);
     x*z/(d1*d2*sqrt(r^2-z^2)), y*z/(d1*d2*sqrt(r^2-z^2)), 1/sqrt(r^2-z^2)*(z^2/(d1*d2)+Upsilon*z*(d2-d1)/(d1*d2*(d1+d2))-1); 
     -y/(x^2+y^2),           x/(x^2+y^2),            0];