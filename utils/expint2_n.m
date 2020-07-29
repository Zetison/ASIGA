function I = expint2_n(omega,n)
% Computes integral(y^n exp(i*omega*y) dy, y=-1..1)

sinw = sin(omega);
cosw = cos(omega);
I = 0;
mfact = 1;
for m = 0:n
    if mod(m,2) % m is odd
        I = I + (-1i*omega).^(m-n-1)/mfact*cosw;
    else
        I = I + 1i*(-1i*omega).^(m-n-1)/mfact*sinw;
    end
    mfact = (m+1)*mfact;
end
I = -2*mfact/(m+1)*I;