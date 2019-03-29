function nurbs = degenerateIGAtoFEM2(nurbs)

coeffs = nurbs.coeffs;

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
        
        for i = 1:p_xi:n_xi-p_xi
            for j = 1:p_eta:n_eta
                for l = 1:p_zeta:n_zeta
                    P_1 = coeffs(1:3,i,j,l);
                    P_end = coeffs(1:3,i+p_xi,j,l);
                    Q = zeros(3,p_xi-1);
                    
                    for ii = 1:p_xi-1
                        Q(:,ii) = (p_xi-ii)*P_1 + ii*P_end;
                        Q(:,ii) = Q(:,ii)/norm(Q(:,ii));
                        Q(:,ii) = Q(:,ii)*norm(P_1);
                    end
                    Q = reshape(Q,3*(p_xi-1),1);
                    C = zeros(p_xi-1, p_xi+1);
                    for m = 1:p_xi-1
                        if i < n_xi
                            i_tilde = i + p_xi - mod(i-1,p_xi);
                        else
                            i_tilde = n_xi;
                        end
                        if j < n_eta
                            j_tilde = j + p_eta - mod(j-1,p_eta);
                        else
                            j_tilde = n_eta;
                        end
                        if l < n_zeta
                            l_tilde = l + p_zeta - mod(l-1,p_zeta);
                        else
                            l_tilde = n_zeta;
                        end
                        xi_tilde = ((p_xi-m)*t_xi(i_tilde) + m*t_xi(i_tilde+1))/p_xi;
                        eta_tilde = t_eta(j_tilde);
                        zeta_tilde = t_zeta(l_tilde);
                        B_xi = Bspline_basisDers(i_tilde, xi_tilde, p_xi, t_xi);
                        B_eta = Bspline_basisDers(j_tilde, eta_tilde, p_eta, t_eta);
                        B_zeta = Bspline_basisDers(l_tilde, zeta_tilde, p_zeta, t_zeta);
                        C(m,:) = B_xi*B_eta(mod(j-1,p_eta)+1)*B_zeta(mod(l-1,p_zeta)+1);
                    end
                    
                    A = kron(C(:,2:end-1), eye(3));
                    P = A\(Q-kron(C(:,1),P_1)-kron(C(:,end),P_end));
                    P = reshape(P,3,p_xi-1);
                    for ii = 1:p_xi-1
                        coeffs(1:3,i+ii,j,l) = P(:,ii);
                    end
                end
            end
        end
        for i = 1:n_xi
            for j = 1:p_eta:n_eta-p_eta
                for l = 1:p_zeta:n_zeta
                    P_1 = coeffs(1:3,i,j,l);
                    P_1 = P_1/norm(P_1)*norm(coeffs(1:3,1,1,l));
                    P_end = coeffs(1:3,i,j+p_eta,l);
                    P_end = P_end/norm(P_end)*norm(coeffs(1:3,1,1,l));
                    Q = zeros(3,p_eta-1);
                    
                    
                    for jj = 1:p_eta-1
                        Q(:,jj) = (p_eta-jj)*P_1 + jj*P_end;
                        Q(:,jj) = Q(:,jj)/norm(Q(:,jj));
                        Q(:,jj) = Q(:,jj)*norm(P_1);
                    end
                    Q = reshape(Q,3*(p_eta-1),1);
                    C = zeros(p_eta-1, p_eta+1);
                    for m = 1:p_eta-1
                        if i < n_xi
                            i_tilde = i + p_xi - mod(i-1,p_xi);
                        else
                            i_tilde = n_xi;
                        end
                        if j < n_eta
                            j_tilde = j + p_eta - mod(j-1,p_eta);
                        else
                            j_tilde = n_eta;
                        end
                        if l < n_zeta
                            l_tilde = l + p_zeta - mod(l-1,p_zeta);
                        else
                            l_tilde = n_zeta;
                        end
                        xi_tilde = t_xi(i_tilde);
                        eta_tilde = ((p_eta-m)*t_eta(j_tilde) + m*t_eta(j_tilde+1))/p_eta;
                        zeta_tilde = t_zeta(l_tilde);
                        B_xi = Bspline_basisDers(i_tilde, xi_tilde, p_xi, t_xi);
                        B_eta = Bspline_basisDers(j_tilde, eta_tilde, p_eta, t_eta);
                        B_zeta = Bspline_basisDers(l_tilde, zeta_tilde, p_zeta, t_zeta);
                        C(m,:) = B_xi(mod(i-1,p_xi)+1)*B_eta*B_zeta(mod(l-1,p_zeta)+1);
                    end
                    
                    A = kron(C(:,2:end-1), eye(3));
                    P = A\(Q-kron(C(:,1),P_1)-kron(C(:,end),P_end));
                    P = reshape(P,3,p_eta-1);
                    for jj = 1:p_eta-1
                        coeffs(1:3,i,j+jj,l) = P(:,jj);
                    end
%                     if i == 2 && j == 3 && l == 1
%                         keyboard
%                     end
                end
            end
        end
        for i = 1:n_xi
            for j = 1:n_eta
                for l = 1:p_zeta:n_zeta-p_zeta
                    for kk = 1:p_zeta-1
                        coeffs(1:3,i,j,l+kk) = ((p_zeta-kk)*coeffs(1:3,i,j,l) + kk*coeffs(1:3,i,j,l+p_zeta))/p_zeta;
                    end
                end
            end
        end   
        
        
        for i = 1:n_xi
            for j = 1:n_eta
                for l = 1:n_zeta
                    coeffs(4,i,j,l) = 1;
                end
            end
        end
end

nurbs.coeffs = coeffs;


