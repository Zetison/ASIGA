function v = evaluateSpline(solid, evaluationPt)

error('Depricated, use evaluateNURBS() instead')

if strcmp(solid.type, '3Dvolume')
    xi = evaluationPt(1);
    eta = evaluationPt(2);
    zeta = evaluationPt(3);
    
    p = solid.degree(1);
    q = solid.degree(2);
    r = solid.degree(3);

    n = solid.number(1);
    m = solid.number(2);
    l = solid.number(3);


    Xi = solid.knots{1};
    Eta = solid.knots{2};
    Zeta = solid.knots{3};

    P = solid.coeffs;

    i1 = findKnotSpan(n, p, xi, Xi);
    i2 = findKnotSpan(m, q, eta, Eta);
    i3 = findKnotSpan(l, r, zeta, Zeta);


    N = Bspline_basis(i1, xi, p, Xi);
    M = Bspline_basis(i2, eta, q, Eta);
    L = Bspline_basis(i3, zeta, r, Zeta);

    v = zeros(3,1);

    for k3 = 1:r+1
        A3 = i3 - r + k3 - 1;
        for k2 = 1:q+1
            A2 = i2 - q + k2 - 1;
            for k1 = 1:p+1
                A1 = i1 - p + k1 - 1;
                v = v + N(k1)*M(k2)*L(k3)*P(1:3, A1, A2, A3);
            end
        end
    end
elseif strcmp(solid.type, '2Dcurve') 
    xi = evaluationPt;
    n = solid.number;
    p = solid.degree;
    Xi = solid.knots;
    P = solid.coeffs;
    i1 = findKnotSpan(n, p, xi, Xi);

    N = Bspline_basis(i1, xi, p, Xi);

    v = zeros(2,1);

    for i = 1:p+1
        A = i1-p+i-1;
        v = v + N(i)*P(1:2,A);
    end
end