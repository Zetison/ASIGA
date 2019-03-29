function nurbs = degenerateIGAtoFEM(nurbs,coreMethod)

P = nurbs.coeffs;

switch nurbs.type
    case '3Dvolume'
        switch coreMethod
            case 'hp_FEM'
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

                i_arr = 1:p_xi:n_xi-1;
                temp = zeros(p_xi,length(i_arr));
                parfor i_p = 1:length(i_arr)        
                    i_pf = i_arr(i_p);

                    dXdxi = @(xi) evaluateNURBS_deriv2(nurbs, [xi, (t_eta(p_eta+1)+t_eta(p_eta+2))/2, (t_zeta(p_zeta+1)+t_zeta(p_zeta+2))/2], 'xi');
                    integrand = @(xi) norm(dXdxi(xi));
                    L_xi = @(xi) quadInt(integrand,t_xi(i_pf+p_xi),xi);
                    dFdxi = @(xi) integrand(xi);

                    temp2 = zeros(p_xi-1,1);
                    for i_tilde = i_pf+1:i_pf+p_xi-1
                        F = @(xi) L_xi(xi) - (i_tilde-i_pf)/p_xi*L_xi(t_xi(i_pf+p_xi+1));
                        temp2(i_tilde-i_pf) = newtonsMethod(F,dFdxi,t_xi(i_pf+p_xi),100,1e-14);
                    end
                    temp(:,i_p) = [t_xi(i_pf+p_xi); temp2];
                end
                xi_tilde(1:end-1) = reshape(temp,n_xi-1,1);
                xi_tilde(end) = t_xi(end);

                j_arr = 1:p_eta:n_eta-1;
                temp = zeros(p_eta,length(j_arr));
                parfor j_p = 1:length(j_arr)    
                    j_pf = j_arr(j_p);   

                    dXdeta = @(eta) evaluateNURBS_deriv2(nurbs, [(t_xi(p_xi+1)+t_xi(p_xi+2))/2, eta, (t_zeta(p_zeta+1)+t_zeta(p_zeta+2))/2], 'eta');
                    integrand = @(eta) norm(dXdeta(eta));
                    L_eta = @(eta) quadInt(integrand,t_eta(j_pf+p_eta),eta);
                    dFdeta = @(eta) integrand(eta);

                    temp2 = zeros(p_eta-1,1);
                    for j_tilde = j_pf+1:j_pf+p_eta-1
                        F = @(eta) L_eta(eta) - (j_tilde-j_pf)/p_eta*L_eta(t_eta(j_pf+p_eta+1));
                        temp2(j_tilde-j_pf) = newtonsMethod(F,dFdeta,t_eta(j_pf+p_eta),100,1e-14);
                    end
                    temp(:,j_p) = [t_eta(j_pf+p_eta); temp2];
                end
                eta_tilde(1:end-1) = reshape(temp,n_eta-1,1);
                eta_tilde(end) = t_eta(end);

                l_arr = 1:p_zeta:n_zeta-1;
                temp = zeros(p_zeta,length(l_arr));
                parfor l_p = 1:length(l_arr)
                    l_pf = l_arr(l_p);   

                    dXdzeta = @(zeta) evaluateNURBS_deriv2(nurbs, [(t_xi(p_xi+1)+t_xi(p_xi+2))/2, (t_eta(p_eta+1)+t_eta(p_eta+2))/2, zeta], 'zeta');
                    integrand = @(zeta) norm(dXdzeta(zeta));
                    L_zeta = @(zeta) quadInt(integrand,t_zeta(l_pf+p_zeta),zeta);
                    dFdzeta = @(zeta) integrand(zeta);

                    temp2 = zeros(p_zeta-1,1);
                    for l_tilde = l_pf+1:l_pf+p_zeta-1
                        F = @(zeta) L_zeta(zeta) - (l_tilde-l_pf)/p_zeta*L_zeta(t_zeta(l_pf+p_zeta+1));
                        temp2(l_tilde-l_pf) = newtonsMethod(F,dFdzeta,t_zeta(l_pf+p_zeta),100,1e-14);
                    end
                    temp(:,l_p) = [t_zeta(l_pf+p_zeta); temp2];
                end
                zeta_tilde(1:end-1) = reshape(temp,n_zeta-1,1);
                zeta_tilde(end) = t_zeta(end);

                parfor i = 1:n_xi
                    for j = 1:n_eta
                        for l = 1:n_zeta
                            if ~(mod(l-1,p_zeta)== 0 && mod(j-1,p_eta) == 0 && mod(i-1,p_xi) == 0)
                                Q(:,i,j,l) = evaluateNURBS(nurbs,[xi_tilde(i),eta_tilde(j),zeta_tilde(l)]);
                            end                                
                        end
                    end
                end


                B_xi_arr = zeros(n_xi,p_xi+1);
                i_xi_arr = zeros(n_xi,1);
                parfor i = 1:n_xi
                    i_xi_p = findKnotSpan(n_xi, p_xi, xi_tilde(i), t_xi)        
                    B_xi_arr(i,:) = Bspline_basisDers(i_xi_p, xi_tilde(i), p_xi, t_xi);
                    i_xi_arr(i) = i_xi_p;
                end

                B_eta_arr = zeros(n_eta,p_eta+1);
                j_eta_arr = zeros(n_eta,1);
                parfor j = 1:n_eta
                    j_eta_p = findKnotSpan(n_eta, p_eta, eta_tilde(j), t_eta)        
                    B_eta_arr(j,:) = Bspline_basisDers(j_eta_p, eta_tilde(j), p_eta, t_eta);
                    j_eta_arr(j) = j_eta_p;
                end

                B_zeta_arr = zeros(n_zeta,p_zeta+1);
                l_zeta_arr = zeros(n_zeta,1);
                parfor l = 1:n_zeta
                    l_zeta_p = findKnotSpan(n_zeta, p_zeta, zeta_tilde(l), t_zeta)        
                    B_zeta_arr(l,:) = Bspline_basisDers(l_zeta_p, zeta_tilde(l), p_zeta, t_zeta);
                    l_zeta_arr(l) = l_zeta_p;
                end

                index = zeros(n_xi*n_eta*n_zeta,3);
                counter = 1;
                for l = 1:n_zeta
                    for j = 1:n_eta
                        for i = 1:n_xi   
                            index(counter,:) = [i, j, l];
                            counter = counter + 1;
                        end
                    end
                end

                values = zeros(n_xi*n_eta*n_zeta,(p_xi+1)*(p_eta+1)*(p_zeta+1));
                n_idx = zeros(n_xi*n_eta*n_zeta,(p_xi+1)*(p_eta+1)*(p_zeta+1));
                m_idx = zeros(n_xi*n_eta*n_zeta,(p_xi+1)*(p_eta+1)*(p_zeta+1));
                parfor a = 1:n_xi*n_eta*n_zeta
                    i = index(a,1);
                    j = index(a,2);
                    l = index(a,3);
                    l_zeta   = l_zeta_arr(l);
                    B_zeta = B_zeta_arr(l,:);
                    j_eta   = j_eta_arr(j);
                    B_eta = B_eta_arr(j,:);
                    i_xi   = i_xi_arr(i);
                    B_xi = B_xi_arr(i,:);
                    m = i + (j-1)*n_xi + (l-1)*n_xi*n_eta;

                    temp = zeros(1,(p_xi+1)*(p_eta+1)*(p_zeta+1));
                    temp2 = zeros(1,(p_xi+1)*(p_eta+1)*(p_zeta+1));
                    counter = 1;
                    for ll = 1:p_zeta+1
                        l_tilde = l_zeta - p_zeta + ll - 1;
                        for jj = 1:p_eta+1
                            j_tilde = j_eta - p_eta + jj - 1;
                            for ii = 1:p_xi+1
                                i_tilde = i_xi - p_xi + ii - 1;
                                n = i_tilde + (j_tilde-1)*n_xi + (l_tilde-1)*n_xi*n_eta;

                                temp(counter) = B_xi(ii)*B_eta(jj)*B_zeta(ll);
                                temp2(counter) = n;
                                counter = counter + 1;
                            end
                        end
                    end
                    values(a,:) = temp;
                    n_idx(a,:) = temp2;
                    m_idx(a,:) = ones(1,(p_xi+1)*(p_eta+1)*(p_zeta+1))*m
                end
                A = sparse(m_idx,n_idx,values);

                Q = reshape(Q,3,n_xi*n_eta*n_zeta);
                P_tilde = zeros(size(Q));
                P_tilde(1,:) = (A\Q(1,:)')';
                P_tilde(2,:) = (A\Q(2,:)')';
                P_tilde(3,:) = (A\Q(3,:)')';
                P(1:3,:,:,:) = reshape(P_tilde,3,n_xi,n_eta,n_zeta);
                P(4,:,:,:) = 1;
            case 'h_FEM'
                n_xi = nurbs.number(1);
                n_eta = nurbs.number(2);
                n_zeta = nurbs.number(3);
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                p_zeta = nurbs.degree(3);
                for i = 1:p_xi:n_xi-p_xi
                    for j = 1:p_eta:n_eta
                        for k = 1:p_zeta:n_zeta
                            for ii = 1:p_xi-1
                                P(1:3,i+ii,j,k) = ((p_xi-ii)*P(1:3,i,j,k) + ii*P(1:3,i+p_xi,j,k))/p_xi;
                            end
                        end
                    end
                end
                for i = 1:n_xi
                    for j = 1:p_eta:n_eta-p_eta
                        for k = 1:p_zeta:n_zeta
                            for jj = 1:p_eta-1
                                P(1:3,i,j+jj,k) = ((p_eta-jj)*P(1:3,i,j,k) + jj*P(1:3,i,j+p_eta,k))/p_eta;
                            end
                        end
                    end
                end
                for i = 1:n_xi
                    for j = 1:n_eta
                        for k = 1:p_zeta:n_zeta-p_zeta
                            for kk = 1:p_zeta-1
                                P(1:3,i,j,k+kk) = ((p_zeta-kk)*P(1:3,i,j,k) + kk*P(1:3,i,j,k+p_zeta))/p_zeta;
                            end
                        end
                    end
                end   
                for i = 1:n_xi
                    for j = 1:n_eta
                        for k = 1:n_zeta
                            P(4,i,j,k) = 1;
                        end
                    end
                end
            case 'linear_FEM'
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                p_zeta = nurbs.degree(3);
                
                P = P(:,1:p_xi:end,1:p_eta:end,1:p_zeta:end);
                P(4,:,:,:) = 1;
                nurbs.knots{1} = [0, unique(nurbs.knots{1}), 1];
                nurbs.knots{2} = [0, unique(nurbs.knots{2}), 1];
                nurbs.knots{3} = [0, unique(nurbs.knots{3}), 1];
                nurbs.degree = [1,1,1];
                nurbs.number = [length(nurbs.knots{1})-2,...
                                length(nurbs.knots{2})-2,...
                                length(nurbs.knots{3})-2];
        end
    case '3Dsurface'
        switch coreMethod
            case 'hp_FEM'
                n_xi = nurbs.number(1);
                n_eta = nurbs.number(2);
                
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                
                t_xi = nurbs.knots{1};
                t_eta = nurbs.knots{2};
                
                Q = zeros(size(P)-[1 0 0]);
                xi_tilde = zeros(n_xi,1);
                eta_tilde = zeros(n_eta,1);
                
                for i = 1:p_xi:n_xi
                    for j = 1:p_eta:n_eta
                        Q(:,i,j) = P(1:3,i,j);
                    end
                end

                i_arr = 1:p_xi:n_xi-1;
                temp = zeros(p_xi,length(i_arr));
                parfor i_p = 1:length(i_arr)
                    i_pf = i_arr(i_p);

                    dXdxi = @(xi) evaluateNURBS_deriv2(nurbs, [xi, (t_eta(p_eta+1)+t_eta(p_eta+2))/2], 'xi');
                    integrand = @(xi) norm(dXdxi(xi));
                    L_xi = @(xi) quadInt(integrand,t_xi(i_pf+p_xi),xi);
                    dFdxi = @(xi) integrand(xi);

                    temp2 = zeros(p_xi-1,1);
                    for i_tilde = i_pf+1:i_pf+p_xi-1
                        F = @(xi) L_xi(xi) - (i_tilde-i_pf)/p_xi*L_xi(t_xi(i_pf+p_xi+1));
                        temp2(i_tilde-i_pf) = newtonsMethod(F,dFdxi,t_xi(i_pf+p_xi),100,1e-14);
                    end
                    temp(:,i_p) = [t_xi(i_pf+p_xi); temp2];
                end
                xi_tilde(1:end-1) = reshape(temp,n_xi-1,1);
                xi_tilde(end) = t_xi(end);

                j_arr = 1:p_eta:n_eta-1;
                temp = zeros(p_eta,length(j_arr));
                parfor j_p = 1:length(j_arr)
                    j_pf = j_arr(j_p);   

                    dXdeta = @(eta) evaluateNURBS_deriv2(nurbs, [(t_xi(p_xi+1)+t_xi(p_xi+2))/2, eta], 'eta');
                    integrand = @(eta) norm(dXdeta(eta));
                    L_eta = @(eta) quadInt(integrand,t_eta(j_pf+p_eta),eta);
                    dFdeta = @(eta) integrand(eta);

                    temp2 = zeros(p_eta-1,1);
                    for j_tilde = j_pf+1:j_pf+p_eta-1
                        F = @(eta) L_eta(eta) - (j_tilde-j_pf)/p_eta*L_eta(t_eta(j_pf+p_eta+1));
                        temp2(j_tilde-j_pf) = newtonsMethod(F,dFdeta,t_eta(j_pf+p_eta),100,1e-14);
                    end
                    temp(:,j_p) = [t_eta(j_pf+p_eta); temp2];
                end
                eta_tilde(1:end-1) = reshape(temp,n_eta-1,1);
                eta_tilde(end) = t_eta(end);

                parfor i = 1:n_xi
                    for j = 1:n_eta
                        if ~(mod(j-1,p_eta) == 0 && mod(i-1,p_xi) == 0)
                            Q(:,i,j) = evaluateNURBS(nurbs,[xi_tilde(i),eta_tilde(j)]);
                        end  
                    end
                end


                B_xi_arr = zeros(n_xi,p_xi+1);
                i_xi_arr = zeros(n_xi,1);
                parfor i = 1:n_xi
                    i_xi_p = findKnotSpan(n_xi, p_xi, xi_tilde(i), t_xi)        
                    B_xi_arr(i,:) = Bspline_basisDers(i_xi_p, xi_tilde(i), p_xi, t_xi);
                    i_xi_arr(i) = i_xi_p;
                end

                B_eta_arr = zeros(n_eta,p_eta+1);
                j_eta_arr = zeros(n_eta,1);
                parfor j = 1:n_eta
                    j_eta_p = findKnotSpan(n_eta, p_eta, eta_tilde(j), t_eta)        
                    B_eta_arr(j,:) = Bspline_basisDers(j_eta_p, eta_tilde(j), p_eta, t_eta);
                    j_eta_arr(j) = j_eta_p;
                end
                
                index = zeros(n_xi*n_eta,2);
                counter = 1;
                for j = 1:n_eta
                    for i = 1:n_xi   
                        index(counter,:) = [i, j];
                        counter = counter + 1;
                    end
                end

                values = zeros(n_xi*n_eta,(p_xi+1)*(p_eta+1));
                n_idx = zeros(n_xi*n_eta,(p_xi+1)*(p_eta+1));
                m_idx = zeros(n_xi*n_eta,(p_xi+1)*(p_eta+1));
                parfor a = 1:n_xi*n_eta
                    i = index(a,1);
                    j = index(a,2);
                    
                    i_xi   = i_xi_arr(i);
                    B_xi = B_xi_arr(i,:);
                    j_eta   = j_eta_arr(j);
                    B_eta = B_eta_arr(j,:);
                    m = i + (j-1)*n_xi;

                    temp = zeros(1,(p_xi+1)*(p_eta+1));
                    temp2 = zeros(1,(p_xi+1)*(p_eta+1));
                    counter = 1;
                    for jj = 1:p_eta+1
                        j_tilde = j_eta - p_eta + jj - 1;
                        for ii = 1:p_xi+1
                            i_tilde = i_xi - p_xi + ii - 1;
                            n = i_tilde + (j_tilde-1)*n_xi;

                            temp(counter) = B_xi(ii)*B_eta(jj);
                            temp2(counter) = n;
                            counter = counter + 1;
                        end
                    end
                    values(a,:) = temp;
                    n_idx(a,:) = temp2;
                    m_idx(a,:) = ones(1,(p_xi+1)*(p_eta+1))*m
                end
                A = sparse(m_idx,n_idx,values);

                Q = reshape(Q,3,n_xi*n_eta);
                P_tilde = zeros(size(Q));
                P_tilde(1,:) = (A\Q(1,:)')';
                P_tilde(2,:) = (A\Q(2,:)')';
                P_tilde(3,:) = (A\Q(3,:)')';
                P(1:3,:,:,:) = reshape(P_tilde,3,n_xi,n_eta);
                P(4,:,:,:) = 1;
            case 'h_FEM'
                n_xi = nurbs.number(1);
                n_eta = nurbs.number(2);
                
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                
                for i = 1:p_xi:n_xi-p_xi
                    for j = 1:p_eta:n_eta
                        for ii = 1:p_xi-1
                            P(1:3,i+ii,j) = ((p_xi-ii)*P(1:3,i,j) + ii*P(1:3,i+p_xi,j))/p_xi;
                        end
                    end
                end
                for i = 1:n_xi
                    for j = 1:p_eta:n_eta-p_eta
                        for jj = 1:p_eta-1
                            P(1:3,i,j+jj) = ((p_eta-jj)*P(1:3,i,j) + jj*P(1:3,i,j+p_eta))/p_eta;
                        end
                    end
                end 
                for i = 1:n_xi
                    for j = 1:n_eta
                        P(4,i,j) = 1;
                    end
                end
            case 'linear_FEM'
                p_xi = nurbs.degree(1);
                p_eta = nurbs.degree(2);
                P = P(:,1:p_xi:end,1:p_eta:end);
                P(4,:,:) = 1;
                nurbs.knots{1} = [0, unique(nurbs.knots{1}), 1];
                nurbs.knots{2} = [0, unique(nurbs.knots{2}), 1];
                nurbs.degree = [1,1];
                nurbs.number = [length(nurbs.knots{1})-2,...
                                length(nurbs.knots{2})-2];
                
        end
end
nurbs.coeffs = P;


