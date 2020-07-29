function nurbs = degenerateIGAtoFEM3(nurbs)
error('Depricated, use degenerateIGAtoFEM instead')

P = nurbs.coeffs;

switch nurbs.type
    case '3Dvolume'
        n_xi = nurbs.number(1);
        n_eta = nurbs.number(2);
        n_zeta = nurbs.number(3);
        p_xi = nurbs.degree(1);
        p_eta = nurbs.degree(2);
        p_zeta = nurbs.degree(3);
        t_xi = nurbs.knots{1};
        t_eta = nurbs.knots{2};
        t_zeta = nurbs.knots{3};
        Q = zeros(size(P)-[1 0 0 0]);
        xi_tilde = zeros(n_xi,1);
        eta_tilde = zeros(n_eta,1);
        zeta_tilde = zeros(n_zeta,1);
        for i = 1:p_xi:n_xi
            for j = 1:p_eta:n_eta
                for l = 1:p_zeta:n_zeta
                    Q(:,i,j,l) = P(1:3,i,j,l);
                end
            end
        end
        tic
        for i = 1:p_xi:n_xi
            xi_tilde(i) = t_xi(i+p_xi);
            if i < n_xi
                dXdxi = @(xi) evaluateNURBS_deriv2(nurbs, [xi, (t_eta(p_eta+1)+t_eta(p_eta+2))/2, (t_zeta(p_zeta+1)+t_zeta(p_zeta+2))/2], 'xi');
                integrand = @(xi) norm(dXdxi(xi));
                L_xi = @(xi) simpson(integrand,t_xi(i+p_xi),xi,100);
                dFdxi = @(xi) integrand(xi);
                for i_tilde = i+1:i+p_xi-1
                    F = @(xi) L_xi(xi) - (i_tilde-i)/p_xi*L_xi(t_xi(i+p_xi+1));
                    xi_tilde(i_tilde) = newtonsMethod(F,dFdxi,t_xi(i+p_xi),100,1e-15);
                end
            end
        end
        for j = 1:p_eta:n_eta
            eta_tilde(j) = t_eta(j+p_eta);

            if j < n_eta
                dXdeta = @(eta) evaluateNURBS_deriv2(nurbs, [(t_xi(p_xi+1)+t_xi(p_xi+2))/2, eta, (t_zeta(p_zeta+1)+t_zeta(p_zeta+2))/2], 'eta');
                integrand = @(eta) norm(dXdeta(eta));
                L_eta = @(eta) simpson(integrand,t_eta(j+p_eta),eta,100);
                dFdeta = @(eta) integrand(eta);
                for j_tilde = j+1:j+p_eta-1
                    F = @(eta) L_eta(eta) - (j_tilde-j)/p_eta*L_eta(t_eta(j+p_eta+1));
                    eta_tilde(j_tilde) = newtonsMethod(F,dFdeta,t_eta(j+p_eta),100,1e-15);
                end
            end
        end
        for l = 1:p_zeta:n_zeta
            zeta_tilde(l) = t_zeta(l+p_zeta);

            if l < n_zeta
                dXdzeta = @(zeta) evaluateNURBS_deriv2(nurbs, [(t_xi(p_xi+1)+t_xi(p_xi+2))/2, (t_eta(p_eta+1)+t_eta(p_eta+2))/2, zeta], 'zeta');
                integrand = @(zeta) norm(dXdzeta(zeta));
                L_zeta = @(zeta) simpson(integrand,t_zeta(l+p_zeta),zeta,100);
                dFdzeta = @(zeta) integrand(zeta);
                for l_tilde = l+1:l+p_zeta-1
                    F = @(zeta) L_zeta(zeta) - (l_tilde-l)/p_zeta*L_zeta(t_zeta(l+p_zeta+1));
                    zeta_tilde(l_tilde) = newtonsMethod(F,dFdzeta,t_zeta(l+p_zeta),100,1e-15);
                end
            end
        end
        toc
        
        for i = 1:n_xi
            for j = 1:n_eta
                for l = 1:n_zeta
                    if ~(mod(l-1,p_zeta)== 0 && mod(j-1,p_eta) == 0 && mod(i-1,p_xi) == 0)
                        Q(:,i,j,l) = evaluateNURBS(nurbs,[xi_tilde(i),eta_tilde(j),zeta_tilde(l)]);
                    end                                
                end
            end
        end
        
        A = sparse(n_xi*n_eta*n_zeta, n_xi*n_eta*n_zeta);
        for l = 1:n_zeta
            for j = 1:n_eta
                for i = 1:n_xi
                    m = i + (j-1)*n_xi + (l-1)*n_xi*n_eta;
                    i_xi   = findKnotSpan(n_xi,   p_xi,   xi_tilde(i),   t_xi);
                    j_eta  = findKnotSpan(n_eta,  p_eta,  eta_tilde(j),  t_eta);
                    l_zeta = findKnotSpan(n_zeta, p_zeta, zeta_tilde(l), t_zeta);
                    
                    B_xi   = Bspline_basisDers(i_xi,   xi_tilde(i),   p_xi,   t_xi);
                    B_eta  = Bspline_basisDers(j_eta,  eta_tilde(j),  p_eta,  t_eta);
                    B_zeta = Bspline_basisDers(l_zeta, zeta_tilde(l), p_zeta, t_zeta);
                    
                    for ll = 1:p_zeta+1
                        l_tilde = l_zeta - p_zeta + ll - 1;
                        for jj = 1:p_eta+1
                            j_tilde = j_eta - p_eta + jj - 1;
                            for ii = 1:p_xi+1
                                i_tilde = i_xi - p_xi + ii - 1;
                                n = i_tilde + (j_tilde-1)*n_xi + (l_tilde-1)*n_xi*n_eta;
                                
                                A(m,n) = B_xi(ii)*B_eta(jj)*B_zeta(ll);
                            end
                        end
                    end
                end
            end
        end
        Q = reshape(Q,3,n_xi*n_eta*n_zeta);
        P_tilde = zeros(size(Q));
        P_tilde(1,:) = (A\Q(1,:)')';
        P_tilde(2,:) = (A\Q(2,:)')';
        P_tilde(3,:) = (A\Q(3,:)')';
        P(1:3,:,:,:) = reshape(P_tilde,3,n_xi,n_eta,n_zeta);
        
        P(4,:,:,:) = 1;
        
end
nurbs.coeffs = P;


