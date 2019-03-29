function [p,dp] = rigidSphereScattering(r,theta,k,R,P_inc,N)

z_R = R*k;
z_r = r*k;
p = zeros(numel(r),numel(k));
dp = zeros(numel(r),numel(k));
Z.R{1,1} = zeros(1,numel(k));
Z.R{2,1} = zeros(1,numel(k));
Z.R{1,2} = bessel_s(0,z_R,1);
Z.R{2,2} = bessel_s(0,z_R,2);
Z.r{1,1} = zeros(numel(r),numel(k));
Z.r{2,1} = zeros(numel(r),numel(k));
Z.r{1,2} = bessel_s(0,z_r,1);
Z.r{2,2} = bessel_s(0,z_r,2);
for n = 0:N
    Z.R{1,1} = Z.R{1,2};
    Z.R{2,1} = Z.R{2,2};
    Z.R{1,2} = bessel_s(n+1,z_R,1);
    Z.R{2,2} = bessel_s(n+1,z_R,2);
    Z.r{1,1} = Z.r{1,2};
    Z.r{2,1} = Z.r{2,2};
    Z.r{1,2} = bessel_s(n+1,z_r,1);
    Z.r{2,2} = bessel_s(n+1,z_r,2);
    
    djR = dbessel_s(n,z_R,1,Z.R);
    dyR = dbessel_s(n,z_R,2,Z.R);
    djR2 = d2bessel_s(n,z_R,1,Z.R);
    dyR2 = d2bessel_s(n,z_R,2,Z.R);
    dhR = djR+1i*dyR;
    dhR2 = djR2+1i*dyR2;
%     Z.r{1,2} = bessel_s(n,z_r,1);
%     Z.r{2,2} = bessel_s(n,z_r,2);
    hr = Z.r{1,1} + 1i*Z.r{2,1};
    djr = dbessel_s(n,z_r,1,Z.r);
    dyr = dbessel_s(n,z_r,2,Z.r);
    dhr = djr + 1i*dyr;
    p = p + 1i^n*(2*n+1)*djR./dhR*legendre_(n,cos(theta)).*hr;
    dp = dp + 1i^n*(2*n+1)*(R*djR2./dhR.*hr - R*djR./dhR.^2.*dhR2.*hr + r*djR./dhR.*dhr)*legendre_(n,cos(theta));
    
end
p = -P_inc*p;
dp = -P_inc*dp;