function Phi = dPhi_kdny(xmy,r,ny,k)
Phi =  Phi_k(r,k)./r.^2.*(1 - 1i*k*r).*sum(xmy.*ny,2);

