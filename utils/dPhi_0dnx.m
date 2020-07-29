function Phi = dPhi_0dnx(xmy,r,nx)
Phi = -Phi_0(r)./r.^2.*(xmy*nx);

