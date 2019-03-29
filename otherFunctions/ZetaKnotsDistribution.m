function D = ZetaKnotsDistribution(coeffs, Zeta_u, zeta_vec)

n = length(coeffs);
p = length(Zeta_u) - n - 1;
D = zeros(size(zeta_vec));

for i = 1:length(zeta_vec)
    zeta_u = zeta_vec(i);
    i1 = findKnotSpan(n, p, zeta_u, Zeta_u);

    N = Bspline_basis(i1, zeta_u, p, Zeta_u, 0);
    
    for k1 = 1:p+1
        A1 = i1 - p + k1 - 1;
        D(i) = D(i) + N(k1)*coeffs(A1);
    end
end