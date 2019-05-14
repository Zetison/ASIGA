function Phi = d2Phi_kdnxdny(xmy,r,nx,ny,k)
Phi = Phi_k(r,k)./r.^2.*((ny*nx).*(1-1i*k*r)+(k^2+3./r.^2.*(1i*k*r-1)).*(xmy*nx).*sum(xmy.*ny,2));