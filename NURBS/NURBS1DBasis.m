function [R, dRdxi] = NURBS1DBasis(xi, p, Xi, weights)

n = length(Xi) - (p+1);

i1 = findKnotSpan(n, p, xi, Xi);

[N, dNdxi] = Bspline_basisDers(i1, xi, p, Xi);

R = zeros(1, p+1);
dRdxi = zeros(1, p+1);

W = 0;
dWdxi = 0;
for k1 = 1:p+1
    weight = weights(k1);

    W       = W       + N(k1)    *weight;
    dWdxi   = dWdxi   + dNdxi(k1)*weight;
end

for k1 = 1:p+1 
    fact = weights(k1)/(W*W);
    
    R(k1) = N(k1)*fact*W;

    dRdxi(k1)   = (dNdxi(k1)*W - N(k1)*dWdxi)*fact;
end

