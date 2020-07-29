function [R, dRdxi, dRdeta, d2Rdxi2, d2Rdeta2] = NURBS2DBasis2(xi, eta, p, q, Xi, Eta, weights)
error('Depricated. Use NURBSbasis instead')

n = length(Xi) - (p+1);
m = length(Eta) - (q+1);

i1 = findKnotSpan(n, p, xi, Xi);
i2 = findKnotSpan(m, q, eta, Eta);

[N, dNdxi, d2Ndxi2] = Bspline_basisDers2(i1, xi, p, Xi);
[M, dMdeta, d2Mdeta2] = Bspline_basisDers2(i2, eta, q, Eta);

R = zeros(1, (p+1)*(q+1));
dRdxi = zeros(1, (p+1)*(q+1));
dRdeta = zeros(1, (p+1)*(q+1));
d2Rdxi2 = zeros(1, (p+1)*(q+1));
d2Rdeta2 = zeros(1, (p+1)*(q+1));

W = 0;
dWdxi = 0;
dWdeta = 0;
d2Wdxi2 = 0;
d2Wdeta2 = 0;

counter = 1;
for k2 = 1:q+1
    for k1 = 1:p+1  
        weight = weights(counter);

        W       = W       + N(k1)    *M(k2)     *weight;
        dWdxi   = dWdxi   + dNdxi(k1)*M(k2)     *weight;
        dWdeta  = dWdeta  + N(k1)    *dMdeta(k2)*weight;
        d2Wdxi2   = d2Wdxi2   + d2Ndxi2(k1)*M(k2)     *weight;
        d2Wdeta2  = d2Wdeta2  + N(k1)    *d2Mdeta2(k2)*weight;
        counter = counter + 1;
    end
end

counter = 1;
for k2 = 1:q+1
    for k1 = 1:p+1    
        fact = weights(counter)/(W*W);

        NM = N(k1)*M(k2);
        R(counter) = NM*fact*W;

        dRdxi(counter)   = (dNdxi(k1)  *M(k2)*W - NM*dWdxi)*fact;
        d2Rdxi2(counter) = ((d2Ndxi2(k1)*W - N(k1)*d2Wdxi2)*fact - 2*fact/W*(dNdxi(k1)*W - N(k1)*dWdxi)*dWdxi)*M(k2);
        dRdeta(counter)  = (dMdeta(k2) *N(k1)*W - NM*dWdeta)*fact;
        d2Rdeta2(counter) = ((d2Mdeta2(k2)*W - M(k2)*d2Wdeta2)*fact - 2*fact/W*(dMdeta(k2)*W - M(k2)*dWdeta)*dWdeta)*N(k1);
        counter = counter + 1;
    end
end

