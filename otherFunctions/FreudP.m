function [P, S] = FreudP(n,r,x)

[B, S] = getFreudCoeffs(n, r);

P = zeros(size(x));
if isa(r,'sym')
    P = vpa(P);
end

for m = 0:n
    P = P + B(m+1)*x.^m;
end
    