function [v, dvdxi, dvdeta, dvdzeta] = evaluateNURBS_deriv(nurbs, parm_pt)
error('Depricated, use evaluateNURBS() instead')

d = size(nurbs.coeffs,1) - 1;

switch nurbs.type
    case '3Dvolume'
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


        [N, dNdxi] = Bspline_basisDers(i1, parm_pt(1), p, Xi);
        [M, dMdeta] = Bspline_basisDers(i2, parm_pt(2), q, Eta);
        [L, dLdzeta] = Bspline_basisDers(i3, parm_pt(3), r, Zeta);

        v = zeros(3,1);
        dvdxi = zeros(3,1);
        dvdeta = zeros(3,1);
        dvdzeta = zeros(3,1);

        W = 0;
        dWdxi = 0;
        dWdeta = 0;
        dWdzeta = 0;

        for k3 = 1:r+1
            A3 = i3 - r + k3 - 1;
            for k2 = 1:q+1
                A2 = i2 - q + k2 - 1;
                for k1 = 1:p+1
                    A1 = i1 - p + k1 - 1;
                    weight = P(4, A1, A2, A3);
                    W       = W       + N(k1)    *M(k2)     *L(k3)      *weight;
                    dWdxi   = dWdxi   + dNdxi(k1)*M(k2)     *L(k3)      *weight;
                    dWdeta  = dWdeta  + N(k1)    *dMdeta(k2)*L(k3)      *weight;
                    dWdzeta = dWdzeta + N(k1)    *M(k2)     *dLdzeta(k3)*weight;
                end
            end
        end

        for k3 = 1:r+1
            A3 = i3 - r + k3 - 1;
            for k2 = 1:q+1
                A2 = i2 - q + k2 - 1;
                for k1 = 1:p+1
                    A1 = i1 - p + k1 - 1;

                    weight = P(4, A1, A2, A3);
                    point = P(1:3, A1, A2, A3);
                    fact = weight/(W*W);

                    NML = N(k1)*M(k2)*L(k3);

                    v = v + N(k1)*M(k2)*L(k3)*point*weight;
                    dvdxi   = dvdxi + (dNdxi(k1)  *M(k2)*L(k3)*W - NML*dWdxi)*fact*point;
                    dvdeta  = dvdeta + (dMdeta(k2) *N(k1)*L(k3)*W - NML*dWdeta)*fact*point;
                    dvdzeta = dvdzeta + (dLdzeta(k3)*N(k1)*M(k2)*W - NML*dWdzeta)*fact*point;
                end
            end
        end
        v = v/W;
    case {'3Dsurface', '2Dsurface'}
        p = nurbs.degree(1);
        q = nurbs.degree(2);

        n = nurbs.number(1);
        m = nurbs.number(2);

        Xi = nurbs.knots{1};
        Eta = nurbs.knots{2};

        P = nurbs.coeffs;

        i1 = findKnotSpan(n, p, parm_pt(1), Xi);
        i2 = findKnotSpan(m, q, parm_pt(2), Eta);


        [N, dNdxi] = Bspline_basisDers(i1, parm_pt(1), p, Xi);
        [M, dMdeta] = Bspline_basisDers(i2, parm_pt(2), q, Eta);

        v = zeros(d,1);
        dvdxi = zeros(d,1);
        dvdeta = zeros(d,1);

        W = 0;
        dWdxi = 0;
        dWdeta = 0;

        for k2 = 1:q+1
            A2 = i2 - q + k2 - 1;
            for k1 = 1:p+1
                A1 = i1 - p + k1 - 1;
                weight = P(d+1, A1, A2);
                W       = W       + N(k1)    *M(k2)     *weight;
                dWdxi   = dWdxi   + dNdxi(k1)*M(k2)     *weight;
                dWdeta  = dWdeta  + N(k1)    *dMdeta(k2)*weight;
            end
        end

        for k2 = 1:q+1
            A2 = i2 - q + k2 - 1;
            for k1 = 1:p+1
                A1 = i1 - p + k1 - 1;

                weight = P(d+1, A1, A2);
                point = P(1:d, A1, A2);
                fact = weight/(W*W);

                NML = N(k1)*M(k2);

                v = v + N(k1)*M(k2)*point*weight;
                dvdxi   = dvdxi + (dNdxi(k1)  *M(k2)*W - NML*dWdxi)*fact*point;
                dvdeta  = dvdeta + (dMdeta(k2) *N(k1)*W - NML*dWdeta)*fact*point;
            end
        end
        v = v/W;
    case {'1Dnurbs', '2Dcurve', '3Dcurve'}
        xi = parm_pt;
        p = nurbs.degree;
        n = nurbs.number;

        Xi = nurbs.knots;
        P = nurbs.coeffs;

        i1 = findKnotSpan(n, p, xi, Xi);

        [N, dNdxi] = Bspline_basisDers(i1, xi, p, Xi);

        v = zeros(d,1);
        dvdxi = zeros(d,1);

        W = 0;
        dWdxi = 0;
        for i = 1:p+1
            A = i1-p+i-1;
            v = v + N(i)*P(1:d,A)*P(end,A);
            W = W + N(i)*P(end,A);
            dWdxi = dWdxi + dNdxi(i)*P(end,A);
        end
        v = v/W;
        for i = 1:p+1
            A = i1-p+i-1;
            weight = P(end, A);
            point = P(1:d, A);
            fact = weight/(W*W);
            dvdxi = dvdxi + (dNdxi(i)*W - N(i)*dWdxi)*fact*point;
        end
end

