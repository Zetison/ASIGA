function Z = sphBessel(n,z,i)


Z = sqrt(pi/2)./sqrt(z).*cylBessel(n+0.5,z,i);