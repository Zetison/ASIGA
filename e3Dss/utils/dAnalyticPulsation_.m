function dp = dAnalyticPulsation_(v,C_n,y,k,n)

Phi_k = @(r) exp(1i*k(1)*r)./(4*pi*r);
dPhi_kdnx = @(xmy,r,ny) Phi_k(r)./r.^2.*(1i*k(1)*r - 1).*sum(xmy.*ny,2);
dp = zeros(size(v,1),1);
for i = 1:size(y,1)
    xms = @(v) [v(:,1)-y(i,1),v(:,2)-y(i,2),v(:,3)-y(i,3)];
    R = @(v) norm2(xms(v));
    dp = dp + C_n(i)*dPhi_kdnx(xms(v),R(v),n);
end