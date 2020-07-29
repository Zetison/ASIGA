function Phi = d2Phi_0dnxdny(xmy,r,nx,ny)
Phi = Phi_0(r)./r.^2.*((ny*nx) - 3./r.^2.*(xmy*nx).*sum(xmy.*ny,2));

