function dZ = dSphBessel(n,z,i)

dZ = n./z.*sphBessel(n,z,i) - sphBessel(n+1,z,i);