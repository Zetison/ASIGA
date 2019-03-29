function h = hankel_s(n,z,type)
%Returns the n'th spherical bessel function of kind "type" evaluated at
%every element in z

if type == 1
    h = bessel_s(n,z,1) + 1i*bessel_s(n,z,2);
else
    h = bessel_s(n,z,1) - 1i*bessel_s(n,z,2);
end