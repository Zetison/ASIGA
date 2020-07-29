function [R, dRdxi, d2Rdxi2, d3Rdxi3] = NURBS1DBasis2(xi, p, Xi, weights)
error('Depricated. Use NURBSbasis instead')

n = length(Xi) - (p+1);

i1 = findKnotSpan(n, p, xi, Xi);

[N, dNdxi, d2Ndxi2, d3Ndxi3] = Bspline_basisDers2(i1, xi, p, Xi);

R = zeros(1, p+1);
dRdxi = zeros(1, p+1);
d2Rdxi2 = zeros(1, p+1);
d3Rdxi3 = zeros(1, p+1);

W = 0;
dWdxi = 0;
d2Wdxi2 = 0;
d3Wdxi3 = 0;

counter = 1;
for k1 = 1:p+1
    weight = weights(counter);

    W       = W       + N(k1)    *weight;
    dWdxi   = dWdxi   + dNdxi(k1)*weight;
    d2Wdxi2   = d2Wdxi2   + d2Ndxi2(k1)*weight;
    d3Wdxi3   = d3Wdxi3   + d3Ndxi3(k1)*weight;
    counter = counter + 1;
end

counter = 1;
for k1 = 1:p+1 
    fact = weights(counter)/(W*W);
    
    R(counter) = N(k1)*fact*W;

    dRdxi(counter)   = (dNdxi(k1)*W - N(k1)*dWdxi)*fact;
    d2Rdxi2(counter) = (d2Ndxi2(k1)*W - N(k1)*d2Wdxi2)*fact - 2*fact/W*(dNdxi(k1)*W - N(k1)*dWdxi)*dWdxi;
    d3Rdxi3(counter) = (d3Ndxi3(k1)*W^2 - 3*d2Ndxi2(k1)*W^2*dWdxi + dNdxi(k1)*(-3*d2Wdxi2*W^2+6*W*dWdxi^2) ...
                        + N(k1)*(-d3Wdxi3*W^2 + 6*d2Wdxi2*dWdxi*W - 6*dWdxi^3))/W^4;
    counter = counter + 1;
end

