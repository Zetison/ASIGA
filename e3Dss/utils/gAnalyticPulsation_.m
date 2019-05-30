function dp = gAnalyticPulsation_(v,C_n,y,k)

Phi_k = @(r) exp(1i*k(1)*r)./(4*pi*r);
dp = zeros(size(v,1),1);
for i = 1:size(y,1)
    xms = @(v) [v(:,1)-y(i,1),v(:,2)-y(i,2),v(:,3)-y(i,3)];
    R = @(v) norm2(xms(v));
    p_i = C_n(i)*Phi_k(R(v));
    dp = dp + elementProd(p_i.*(1i*k - 1./R(v))./R(v), xms(v));
end