function Phi = dPhi_0dny(xmy,r,ny)
Phi = Phi_0(r)./r.^2.*sum(xmy.*ny,2);

