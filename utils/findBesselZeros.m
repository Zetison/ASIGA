function findBesselZeros()


noFuncs = 100;
noZeros = 100;
besselJZeros = zeros(noFuncs,noZeros);
dbesselJZeros = zeros(noFuncs,noZeros);
besselYZeros = zeros(noFuncs,noZeros);
dbesselYZeros = zeros(noFuncs,noZeros);
for n = 1:noFuncs
    I = [0.1,10*noZeros];
    N = 10*noZeros;
    
    besselJZero = findZeros(@(x) besselj(n-1,x),I,N);
    besselJZeros(n,:) = besselJZero(1:noZeros);
    dbesselJZero = findZeros(@(x) (n-1)./x.*besselj(n-1,x) - besselj(n,x),I,N);
    dbesselJZeros(n,:) = dbesselJZero(1:noZeros);
    
    besselYZero = findZeros(@(x) bessely(n-1,x),I,N);
    besselYZeros(n,:) = besselYZero(1:noZeros);
    dbesselYZero = findZeros(@(x) (n-1)./x.*bessely(n-1,x) - bessely(n,x),I,N);
    dbesselYZeros(n,:) = dbesselYZero(1:noZeros);
end
% Compare with https://mathworld.wolfram.com/BesselFunctionZeros.html
besselJZeros(1:6,1:5).'
dbesselJZeros(1:6,1:5).'
% To add trivial solution uncomment the following lines
% besselJZeros = [besselJZeros(1,:); [zeros(noFuncs-1,1), besselJZeros(2:end,1:end-1)]];
% dbesselJZeros = [0, dbesselJZeros(1,1:end-1); dbesselJZeros(2,:); [zeros(noFuncs-2,1), dbesselJZeros(3:end,1:end-1)]];
save('miscellaneous/besselZeros/besselJZeros.mat','besselJZeros')
save('miscellaneous/besselZeros/dbesselJZeros.mat','dbesselJZeros')
save('miscellaneous/besselZeros/besselYZeros.mat','besselYZeros')
save('miscellaneous/besselZeros/dbesselYZeros.mat','dbesselYZeros')


besseljZeros = zeros(noFuncs,noZeros);
dbesseljZeros = zeros(noFuncs,noZeros);
besselyZeros = zeros(noFuncs,noZeros);
dbesselyZeros = zeros(noFuncs,noZeros);
for n = 1:noFuncs
    I = [0.1,10*noZeros];
    N = 10*noZeros;
    
    besseljZero = findZeros(@(x) besselj(n-1+1/2,x).*sqrt(pi./(2*x)),I,N);
    besseljZeros(n,:) = besseljZero(1:noZeros);
    dbesseljZero = findZeros(@(x) ((n-1)./x.*besselj(n-1+1/2,x) - besselj(n+1/2,x)).*sqrt(pi./(2*x)),I,N);
    dbesseljZeros(n,:) = dbesseljZero(1:noZeros);
    
    besselyZero = findZeros(@(x) bessely(n-1+1/2,x).*sqrt(pi./(2*x)),I,N);
    besselyZeros(n,:) = besselyZero(1:noZeros);
    dbesselyZero = findZeros(@(x) ((n-1)./x.*bessely(n-1+1/2,x) - bessely(n+1/2,x)).*sqrt(pi./(2*x)),I,N);
    dbesselyZeros(n,:) = dbesselyZero(1:noZeros);
end
% To add trivial solution uncomment the following lines
% besseljZeros = [besseljZeros(1,:); [zeros(noFuncs-1,1), besseljZeros(2:end,1:end-1)]];
% dbesseljZeros = [0, dbesseljZeros(1,1:end-1); dbesseljZeros(2,:); [zeros(noFuncs-2,1), dbesseljZeros(3:end,1:end-1)]];
save('miscellaneous/besselZeros/besseljZeros.mat','besseljZeros')
save('miscellaneous/besselZeros/dbesseljZeros.mat','dbesseljZeros')
save('miscellaneous/besselZeros/besselyZeros.mat','besselyZeros')
save('miscellaneous/besselZeros/dbesselyZeros.mat','dbesselyZeros')