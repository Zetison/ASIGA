function v = evaluateNURBS(nurbs, parm_pt)

d = size(nurbs.coeffs,1) - 1;

if strcmp(nurbs.type, '3Dvolume')
    p = nurbs.degree(1);
    q = nurbs.degree(2);
    r = nurbs.degree(3);

    n = nurbs.number(1);
    m = nurbs.number(2);
    l = nurbs.number(3);

    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    Zeta = nurbs.knots{3};

    P = nurbs.coeffs;

    i1 = findKnotSpan(n, p, parm_pt(1), Xi);
    i2 = findKnotSpan(m, q, parm_pt(2), Eta);
    i3 = findKnotSpan(l, r, parm_pt(3), Zeta);


    N = Bspline_basis(i1, parm_pt(1), p, Xi, 0);
    M = Bspline_basis(i2, parm_pt(2), q, Eta, 0);
    L = Bspline_basis(i3, parm_pt(3), r, Zeta, 0);

    v = zeros(3,1);

    W = 0;
    for k3 = 1:r+1
        A3 = i3 - r + k3 - 1;
        for k2 = 1:q+1
            A2 = i2 - q + k2 - 1;
            for k1 = 1:p+1
                A1 = i1 - p + k1 - 1;
                v = v + N(k1)*M(k2)*L(k3)*P(1:3, A1, A2, A3)*P(4, A1, A2, A3);
                W = W + N(k1)*M(k2)*L(k3)*P(4, A1, A2, A3);
            end
        end
    end
    v = v/W;
elseif strcmp(nurbs.type, '3Dsurface')
    p = nurbs.degree(1);
    q = nurbs.degree(2);

    n = nurbs.number(1);
    m = nurbs.number(2);

    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};

    P = nurbs.coeffs;

    i1 = findKnotSpan(n, p, parm_pt(1), Xi);
    i2 = findKnotSpan(m, q, parm_pt(2), Eta);


    N = Bspline_basis(i1, parm_pt(1), p, Xi, 0);
    M = Bspline_basis(i2, parm_pt(2), q, Eta, 0);

    v = zeros(3,1);

    W = 0;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            v = v + N(k1)*M(k2)*P(1:3, A1, A2)*P(4, A1, A2);
            W = W + N(k1)*M(k2)*P(4, A1, A2);
        end
    end
    v = v/W;
elseif strcmp(nurbs.type, '2Dsurface')
    p = nurbs.degree(1);
    q = nurbs.degree(2);

    n = nurbs.number(1);
    m = nurbs.number(2);

    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};

    P = nurbs.coeffs;

    i1 = findKnotSpan(n, p, parm_pt(1), Xi);
    i2 = findKnotSpan(m, q, parm_pt(2), Eta);


    N = Bspline_basis(i1, parm_pt(1), p, Xi, 0);
    M = Bspline_basis(i2, parm_pt(2), q, Eta, 0);

    v = zeros(2,1);

    W = 0;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            v = v + N(k1)*M(k2)*P(1:2, A1, A2)*P(3, A1, A2);
            W = W + N(k1)*M(k2)*P(3, A1, A2);
        end
    end
    v = v/W;
elseif strcmp(nurbs.type, '2Dcurve') || strcmp(nurbs.type, '3Dcurve')
    xi = parm_pt;
    p = nurbs.degree;
    n = nurbs.number;
    
    Xi = nurbs.knots{1};
    P = nurbs.coeffs;
    
    i1 = findKnotSpan(n, p, xi, Xi);

    N = Bspline_basis(i1, xi, p, Xi, 0);

    v = zeros(d,1);

    W = 0;
    for i = 1:p+1
        A = i1-p+i-1;
        v = v + N(i)*P(1:d,A)*P(d+1,A);
        W = W + N(i)*P(d+1,A);
    end
    v = v/W;    
elseif strcmp(nurbs.type, '1Dnurbs')
    xi = parm_pt;
    p = nurbs.degree;
    n = nurbs.number;
    
    Xi = nurbs.knots{1};
    P = nurbs.coeffs;
    
    i1 = findKnotSpan(n, p, xi, Xi);

    N = Bspline_basis(i1, xi, p, Xi, 0);

    v = 0;

    W = 0;
    for i = 1:p+1
        A = i1-p+i-1;
        v = v + N(i)*P(1,A)*P(2,A);
        W = W + N(i)*P(2,A);
    end
    v = v/W;
end