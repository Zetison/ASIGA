function f = bessel_c(nu,z,type)

if type == 1 % besselj
    f = besselj(nu,z);
else % bessely
    f = bessely(nu,z);
end