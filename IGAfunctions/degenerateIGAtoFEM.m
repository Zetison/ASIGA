function task = degenerateIGAtoFEM(task)
for i_varCol = 1:numel(task.varCol) % assume coreMethod to be the same in all domains
    nurbsPatches = task.varCol{i_varCol}.nurbs;
    if ~iscell(nurbsPatches)
        nurbsPatches = {nurbsPatches};
    end
    if strcmp(task.misc.coreMethod,'SEM')
        nurbsPatches = explodeNURBS(nurbsPatches,1);
        d_p = nurbsPatches{1}.d_p;
        if d_p >= 2
            nurbsPatches = explodeNURBS(nurbsPatches,2);
        end
        if d_p == 3
            nurbsPatches = explodeNURBS(nurbsPatches,3);
        end
    end
    for patch = 1:numel(nurbsPatches)
        nurbs = nurbsPatches{patch};
        P = nurbs.coeffs;

        switch nurbs.d_p
            case 3
                switch task.misc.coreMethod
                    case 'SEM'
                        Nxi = nurbs.number(1);
                        Neta = nurbs.number(2);
                        Nzeta = nurbs.number(3);
                        [tXi,rhoXi] = gaussLobattoLegendreQuad(Nxi);
                        if Neta == Nxi
                            tEta = tXi;
                            rhoEta = rhoXi;
                        else
                            [tEta,rhoEta] = gaussLobattoLegendreQuad(Neta);                        
                        end
                        if Nzeta == Nxi
                            tZeta = tXi;
                            rhoZeta = rhoXi;
                        elseif Nzeta == Neta
                            tZeta = tEta;
                            rhoZeta = rhoEta;
                        else
                            [tZeta,rhoZeta] = gaussLobattoLegendreQuad(Nzeta);                        
                        end
                        P = zeros(4,Nxi,Neta,Nzeta);
                        tXit = parent2ParametricSpace([0,1],tXi);
                        tEtat = parent2ParametricSpace([0,1],tEta);
                        tZetat = parent2ParametricSpace([0,1],tZeta);
                        [xi,eta,zeta] = ndgrid(tXit,tEtat,tZetat);
                        noDofs = Nxi*Neta*Nzeta;
                        temp = evaluateNURBS(nurbs, [reshape(xi,noDofs,1),reshape(eta,noDofs,1),reshape(zeta,noDofs,1)]);
    %                     parfor l = 1:Nzeta
    %                         Ptemp = ones(4,Nxi,Neta);
    %                         for j = 1:Neta
    %                             for i = 1:Nxi
    %                                 Ptemp(1:3,i,j) = evaluateNURBS(nurbs,[tXit(i),tEtat(j),tZetat(l)]);
    %                             end
    %                         end
    %                         P(:,:,:,l) = Ptemp;
    %                     end
                        P(1:3,:,:,:) = reshape(temp.',3,Nxi,Neta,Nzeta);
                        nurbs.rho = {rhoXi,rhoEta,rhoZeta};
                        nurbs.GLL = {tXi,tEta,tZeta};

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
                        for i = 1:p_xi:n_xi
                            for j = 1:p_eta:n_eta
                                for l = 1:p_zeta:n_zeta
                                    Q(:,i,j,l) = P(1:3,i,j,l);
                                end
                            end
                        end
    % 
    %                     xi_tilde = splinesGLL(t_xi,p_xi+1);
    %                     eta_tilde = splinesGLL(t_eta,p_eta+1);
    %                     zeta_tilde = splinesGLL(t_zeta,p_zeta+1);
                        xi_tilde = aveknt(t_xi, p_xi+1);
                        eta_tilde = aveknt(t_eta, p_eta+1);
                        zeta_tilde = aveknt(t_zeta, p_zeta+1);

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
                            B_xi_arr(i,:) = BsplineBasis(i_xi_p, xi_tilde(i), p_xi, t_xi);
                            i_xi_arr(i) = i_xi_p;
                        end

                        B_eta_arr = zeros(n_eta,p_eta+1);
                        j_eta_arr = zeros(n_eta,1);
                        parfor j = 1:n_eta
                            j_eta_p = findKnotSpan(n_eta, p_eta, eta_tilde(j), t_eta)        
                            B_eta_arr(j,:) = BsplineBasis(j_eta_p, eta_tilde(j), p_eta, t_eta);
                            j_eta_arr(j) = j_eta_p;
                        end

                        B_zeta_arr = zeros(n_zeta,p_zeta+1);
                        l_zeta_arr = zeros(n_zeta,1);
                        parfor l = 1:n_zeta
                            l_zeta_p = findKnotSpan(n_zeta, p_zeta, zeta_tilde(l), t_zeta)        
                            B_zeta_arr(l,:) = BsplineBasis(l_zeta_p, zeta_tilde(l), p_zeta, t_zeta);
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
                            m_idx(a,:) = ones(1,(p_xi+1)*(p_eta+1)*(p_zeta+1))*m;
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
            case 2
                switch task.misc.coreMethod
                    case 'SEM'
                        Nxi = nurbs.number(1);
                        Neta = nurbs.number(2);
                        [tXi,rhoXi] = gaussLobattoLegendreQuad(Nxi);
                        if Neta == Nxi
                            tEta = tXi;
                            rhoEta = rhoXi;
                        else
                            [tEta,rhoEta] = gaussLobattoLegendreQuad(Neta);                        
                        end
                        P = zeros(4,Nxi,Neta);
                        tXit = parent2ParametricSpace([0,1],tXi);
                        tEtat = parent2ParametricSpace([0,1],tEta);
                        [xi,eta] = ndgrid(tXit,tEtat);
                        noDofs = Nxi*Neta;
                        temp = evaluateNURBS_2ndDeriv(nurbs, [reshape(xi,noDofs,1),reshape(eta,noDofs,1)]);
    %                     parfor l = 1:Nzeta
    %                         Ptemp = ones(4,Nxi,Neta);
    %                         for j = 1:Neta
    %                             for i = 1:Nxi
    %                                 Ptemp(1:3,i,j) = evaluateNURBS(nurbs,[tXit(i),tEtat(j),tZetat(l)]);
    %                             end
    %                         end
    %                         P(:,:,:,l) = Ptemp;
    %                     end
                        P(1:3,:,:) = reshape(temp.',3,Nxi,Neta);
                        nurbs.rho = {rhoXi,rhoEta};
                        nurbs.GLL = {tXi,tEta};

                    case 'hp_FEM'
                        n_xi = nurbs.number(1);
                        n_eta = nurbs.number(2);

                        p_xi = nurbs.degree(1);
                        p_eta = nurbs.degree(2);

                        t_xi = nurbs.knots{1};
                        t_eta = nurbs.knots{2};

                        Q = zeros(size(P)-[1 0 0]);

                        for i = 1:p_xi:n_xi
                            for j = 1:p_eta:n_eta
                                Q(:,i,j) = P(1:3,i,j);
                            end
                        end
    % 
    %                     xi_tilde = splinesGLL(t_xi,p_xi+1);
    %                     eta_tilde = splinesGLL(t_eta,p_eta+1);
                        xi_tilde = aveknt(t_xi, p_xi+1);
                        eta_tilde = aveknt(t_eta, p_eta+1);

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
                            B_xi_arr(i,:) = BsplineBasis(i_xi_p, xi_tilde(i), p_xi, t_xi);
                            i_xi_arr(i) = i_xi_p;
                        end

                        B_eta_arr = zeros(n_eta,p_eta+1);
                        j_eta_arr = zeros(n_eta,1);
                        parfor j = 1:n_eta
                            j_eta_p = findKnotSpan(n_eta, p_eta, eta_tilde(j), t_eta)        
                            B_eta_arr(j,:) = BsplineBasis(j_eta_p, eta_tilde(j), p_eta, t_eta);
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
        nurbsPatches{patch} = nurbs;
    end
    task.varCol{i_varCol}.nurbs = nurbsPatches;
end
