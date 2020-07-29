function E = incompleteEllipticIntegral(phi,m)

E = integral(@(theta) sqrt(1-m*sin(theta).^2),0,phi);