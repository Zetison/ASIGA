
G = @(r) exp(1i*k*r)./(4*pi*r);
dGdn = @(xms,r,n) G(r).*(1i*k*r - 1)./r.^2 .*dot3(xms, n);
Gt = @(r) 1./(4*pi*r);
dGtdn = @(xms,r,n) -Gt(r)./r.*dot3(xms, n);