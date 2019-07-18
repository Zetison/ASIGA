function varargout = evaluateNURBS_2ndDeriv(nurbs, parm_pt)

d = size(nurbs.coeffs,1) - 1;
computeDers = nargout > 1;
switch nurbs.type
    case '3Dvolume'
        p_xi = nurbs.degree(1);
        p_eta = nurbs.degree(2);
        p_zeta = nurbs.degree(3);

        n_xi = nurbs.number(1);
        n_eta = nurbs.number(2);
        n_zeta = nurbs.number(3);

        Xi = nurbs.knots{1};
        Eta = nurbs.knots{2};
        Zeta = nurbs.knots{3};

        P = nurbs.coeffs;

        i1 = findKnotSpan(n_xi, p_xi, parm_pt(1,1), Xi);
        i2 = findKnotSpan(n_eta, p_eta, parm_pt(1,2), Eta);
        i3 = findKnotSpan(n_zeta, p_zeta, parm_pt(1,3), Zeta);


        N = Bspline_basisDers2(i1, parm_pt(:,1), p_xi, Xi);
        M = Bspline_basisDers2(i2, parm_pt(:,2), p_eta, Eta);
        L = Bspline_basisDers2(i3, parm_pt(:,3), p_zeta, Zeta);

        noxi = size(parm_pt,1);
        W = zeros(noxi,1,class(parm_pt));

        for k3 = 1:p_zeta+1
            A3 = i3 - p_zeta + k3 - 1;
            for k2 = 1:p_eta+1
                A2 = i2 - p_eta + k2 - 1;
                for k1 = 1:p_xi+1
                    A1 = i1 - p_xi + k1 - 1;
                    weight = P(d+1, A1, A2, A3);

                    W           = W             + N(:,k1)     .*M(:,k2)     .*L(:,k3)     *weight;
                end
            end
        end

        v = zeros(noxi,d,class(parm_pt));
        counter = 1;
        for k3 = 1:p_zeta+1
            A3 = i3 - p_zeta + k3 - 1;
            for k2 = 1:p_eta+1
                A2 = i2 - p_eta + k2 - 1;
                for k1 = 1:p_xi+1
                    A1 = i1 - p_xi + k1 - 1;
                    weight = P(d+1, A1, A2, A3);   
                    point = P(1:d, A1, A2, A3).';   
                    fact = weight./(W.*W);

                    NML = N(:,k1).*M(:,k2).*L(:,k3);
                    v = v + NML.*W.*fact*point;                
                    counter = counter + 1;
                end
            end
        end
        varargout = cell(6,1);
        varargout{1} = v;
    case {'3Dsurface', '2Dsurface'}
        p_xi = nurbs.degree(1);
        p_eta = nurbs.degree(2);

        n_xi = nurbs.number(1);
        n_eta = nurbs.number(2);

        Xi = nurbs.knots{1};
        Eta = nurbs.knots{2};

        P = nurbs.coeffs;

        i1 = findKnotSpan(n_xi, p_xi, parm_pt(1,1), Xi);
        i2 = findKnotSpan(n_eta, p_eta, parm_pt(1,2), Eta);

        if computeDers
            [N, dNdxi, d2Ndxi2] = Bspline_basisDers2(i1, parm_pt(:,1), p_xi, Xi);
            [M, dMdeta, d2Mdeta2] = Bspline_basisDers2(i2, parm_pt(:,2), p_eta, Eta);
        else
            N = Bspline_basisDers2(i1, parm_pt(:,1), p_xi, Xi);
            M = Bspline_basisDers2(i2, parm_pt(:,2), p_eta, Eta);
        end

        noxi = size(parm_pt,1);
        W = zeros(noxi,1,class(parm_pt));
        dWdxi = zeros(noxi,1,class(parm_pt));
        dWdeta = zeros(noxi,1,class(parm_pt));
        d2Wdxi2 = zeros(noxi,1,class(parm_pt));
        d2Wdeta2 = zeros(noxi,1,class(parm_pt));
        d2Wdxideta = zeros(noxi,1,class(parm_pt));

        for k2 = 1:p_eta+1
            A2 = i2 - p_eta + k2 - 1;
            for k1 = 1:p_xi+1
                A1 = i1 - p_xi + k1 - 1;
                weight = P(d+1, A1, A2);

                W           = W             + N(:,k1)     .*M(:,k2)     *weight;
                if computeDers
                    dWdxi       = dWdxi         + dNdxi(:,k1) .*M(:,k2)     *weight;
                    dWdeta      = dWdeta        + N(:,k1)     .*dMdeta(:,k2)*weight;
                    d2Wdxi2     = d2Wdxi2       + d2Ndxi2(:,k1).*M(:,k2)     *weight;
                    d2Wdeta2    = d2Wdeta2      + N(:,k1)     .*d2Mdeta2(:,k2)*weight;
                    d2Wdxideta  = d2Wdxideta	+ dNdxi(:,k1) .*dMdeta(:,k2)*weight;
                end
            end
        end

        v = zeros(noxi,d,class(parm_pt));
        if computeDers
            dvdxi = zeros(noxi,d,class(parm_pt));
            dvdeta = zeros(noxi,d,class(parm_pt));
            d2vdxi2 = zeros(noxi,d,class(parm_pt));
            d2vdeta2 = zeros(noxi,d,class(parm_pt));
            d2vdxideta = zeros(noxi,d,class(parm_pt));
        end
        counter = 1;
        for k2 = 1:p_eta+1
            A2 = i2 - p_eta + k2 - 1;
            for k1 = 1:p_xi+1
                A1 = i1 - p_xi + k1 - 1;
                weight = P(d+1, A1, A2);   
                point = P(1:d, A1, A2).';   
                fact = weight./(W.*W);

                NM = N(:,k1).*M(:,k2);
                v = v + NM.*W.*fact*point;
                if computeDers
                    dvdxi   = dvdxi + (dNdxi(:,k1)  .*M(:,k2).*W - NM.*dWdxi).*fact*point;
                    dvdeta  = dvdeta + (dMdeta(:,k2) .*N(:,k1).*W - NM.*dWdeta).*fact*point;

                    d2vdxi2  = d2vdxi2  + (d2Ndxi2(:,k1).*W  - N(:,k1).*d2Wdxi2  - 2*(dNdxi(:,k1).*W  - N(:,k1).*dWdxi).*dWdxi./W).*M(:,k2).*fact*point;
                    d2vdeta2 = d2vdeta2 + (d2Mdeta2(:,k2).*W - M(:,k2).*d2Wdeta2 - 2*(dMdeta(:,k2).*W - M(:,k2).*dWdeta).*dWdeta./W).*N(:,k1).*fact*point;
                    d2vdxideta   = d2vdxideta + (-dNdxi(:,k1).*dWdeta.*M(:,k2).*W - N(:,k1).*dWdxi.*dMdeta(:,k2).*W - NM.*d2Wdxideta.*W ...
                                                    + dNdxi(:,k1).*dMdeta(:,k2).*W.^2 + 2*dWdxi.*dWdeta.*NM)./W.*fact*point;
                end
                counter = counter + 1;
            end
        end
        varargout = cell(6,1);
        varargout{1} = v;
        if computeDers
            varargout{2} = dvdxi;
            varargout{3} = dvdeta;
            varargout{4} = d2vdxi2;
            varargout{5} = d2vdxideta;
            varargout{6} = d2vdeta2;
        end
    case {'1Dnurbs', '2Dcurve', '3Dcurve'}
        xi = parm_pt;
        p_xi = nurbs.degree;
        n_xi = nurbs.number;

        Xi = nurbs.knots;
        P = nurbs.coeffs;

        i1 = findKnotSpan(n_xi, p_xi, xi, Xi);

        [N, dNdxi, d2Ndxi2] = Bspline_basisDers2(i1, xi, p_xi, Xi);

        v = zeros(d,1,class(parm_pt));
        dvdxi = zeros(d,1,class(parm_pt));
        
        d2vdxi2 = zeros(d,1,class(parm_pt));

        W = 0;
        dWdxi = 0;
        d2Wdxi2 = 0;
        for i = 1:p_xi+1
            A = i1-p_xi+i-1;
            weight = P(end, A);
            v = v + N(i)*P(1:d,A)*weight;
            W = W + N(i)*P(end,A);
            dWdxi = dWdxi + dNdxi(i)*weight;
            d2Wdxi2 = d2Wdxi2 + d2Ndxi2(i)*weight;
        end
        v = v/W;
        for i = 1:p_xi+1
            A = i1-p_xi+i-1;
            weight = P(end, A);
            point = P(1:d, A);
            fact = weight/(W*W);
            dvdxi = dvdxi + (dNdxi(i)*W - N(i)*dWdxi)*fact*point;
            d2vdxi2  = d2vdxi2  + (d2Ndxi2(i).*W  - N(i).*d2Wdxi2  - 2*(dNdxi(i).*W  - N(i).*dWdxi).*dWdxi./W).*fact*point;
        end
        varargout = cell(6,1);
        varargout{1} = v;
        varargout{2} = dvdxi;
        varargout{3} = d2vdxi2;
end

