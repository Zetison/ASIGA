function Z = cylBessel(n,z,i)

switch i
    case 1
        Z = besselj(n,z);
    case 2
        Z = bessely(n,z);
    case 3
        Z = besselj(n,z) + 1i*bessely(n,z);
    case 4
        Z = besselj(n,z) - 1i*bessely(n,z);
end