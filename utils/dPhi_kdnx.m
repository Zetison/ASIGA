function Phi = dPhi_kdnx(xmy,r,nx,k)
Phi = -Phi_k(r,k)./r.^2.*(1 - 1i*k*r).*(xmy*nx);

