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
        d_p = nurbs.d_p;
        P = nurbs.coeffs;
        degree = nurbs.degree;
        number = nurbs.number;
        t = nurbs.knots;

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
                        number = zeros(1,d_p);
                        t = cell(1,d_p);
                        xi_t = cell(1,d_p);
                        degree = min(nurbs.degree,task.msh.degree);
                        for i = 1:nurbs.d_p
                            t{i} = [0,repelem(unique(nurbs.knots{i}),degree(i)),1];
                            number(i) = numel(t{i})-(degree(i)+1);
                            xi_t{i} = aveknt(t{i}, degree(i)+1);
%                             xi_t{1}{i} = splinesGLL(t{i}, degree(i)+1);
                        end

                        % Evaluate original nurbs at Greville abscissa for interpolation
                        [xi,eta,zeta] = ndgrid(xi_t{1},xi_t{2},xi_t{3});                        
                        Q = reshape(evaluateNURBSvec(nurbs,[xi(:),eta(:),zeta(:)]).',[nurbs.d,number]);


                        B_xi_arr = zeros(number(1),degree(1)+1);
                        i_xi_arr = zeros(number(1),1);
                        parfor i = 1:number(1)
                            i_xi_p = findKnotSpan(number(1), degree(1), xi_t{1}(i), t{1})        
                            B_xi_arr(i,:) = BsplineBasis(i_xi_p, xi_t{1}(i), degree(1), t{1});
                            i_xi_arr(i) = i_xi_p;
                        end

                        B_eta_arr = zeros(number(2),degree(2)+1);
                        j_eta_arr = zeros(number(2),1);
                        parfor j = 1:number(2)
                            j_eta_p = findKnotSpan(number(2), degree(2), xi_t{2}(j), t{2})        
                            B_eta_arr(j,:) = BsplineBasis(j_eta_p, xi_t{2}(j), degree(2), t{2});
                            j_eta_arr(j) = j_eta_p;
                        end

                        B_zeta_arr = zeros(number(3),degree(3)+1);
                        l_zeta_arr = zeros(number(3),1);
                        parfor l = 1:number(3)
                            l_zeta_p = findKnotSpan(number(3), degree(3), xi_t{3}(l), t{3})        
                            B_zeta_arr(l,:) = BsplineBasis(l_zeta_p, xi_t{3}(l), degree(3), t{3});
                            l_zeta_arr(l) = l_zeta_p;
                        end

                        index = zeros(prod(number),3);
                        counter = 1;
                        for l = 1:number(3)
                            for j = 1:number(2)
                                for i = 1:number(1)   
                                    index(counter,:) = [i, j, l];
                                    counter = counter + 1;
                                end
                            end
                        end

                        values = zeros(prod(number),prod(degree+1));
                        n_idx = zeros(prod(number),prod(degree+1));
                        m_idx = zeros(prod(number),prod(degree+1));
                        parfor a = 1:prod(number)
                            i = index(a,1);
                            j = index(a,2);
                            l = index(a,3);
                            l_zeta   = l_zeta_arr(l);
                            B_zeta = B_zeta_arr(l,:);
                            j_eta   = j_eta_arr(j);
                            B_eta = B_eta_arr(j,:);
                            i_xi   = i_xi_arr(i);
                            B_xi = B_xi_arr(i,:);
                            m = i + (j-1)*number(1) + (l-1)*number(1)*number(2);

                            temp = zeros(1,prod(degree+1));
                            temp2 = zeros(1,prod(degree+1));
                            counter = 1;
                            for ll = 1:degree(3)+1
                                l_tilde = l_zeta - degree(3) + ll - 1;
                                for jj = 1:degree(2)+1
                                    j_tilde = j_eta - degree(2) + jj - 1;
                                    for ii = 1:degree(1)+1
                                        i_tilde = i_xi - degree(1) + ii - 1;
                                        n = i_tilde + (j_tilde-1)*number(1) + (l_tilde-1)*number(1)*number(2);

                                        temp(counter) = B_xi(ii)*B_eta(jj)*B_zeta(ll);
                                        temp2(counter) = n;
                                        counter = counter + 1;
                                    end
                                end
                            end
                            values(a,:) = temp;
                            n_idx(a,:) = temp2;
                            m_idx(a,:) = ones(1,prod(degree+1))*m;
                        end
                        A = sparse(m_idx,n_idx,values);

                        Q = reshape(Q,3,prod(number));
                        P = ones([nurbs.d+1,number]);
                        P_tilde = zeros(size(Q));
                        P_tilde(1,:) = (A\Q(1,:)')';
                        P_tilde(2,:) = (A\Q(2,:)')';
                        P_tilde(3,:) = (A\Q(3,:)')';
                        P(1:3,:,:,:) = reshape(P_tilde,3,number(1),number(2),number(3));
                    case 'h_FEM'
                        number(1) = nurbs.number(1);
                        number(2) = nurbs.number(2);
                        number(3) = nurbs.number(3);
                        degree(1) = nurbs.degree(1);
                        degree(2) = nurbs.degree(2);
                        degree(3) = nurbs.degree(3);
                        for i = 1:degree(1):number(1)-degree(1)
                            for j = 1:degree(2):number(2)
                                for k = 1:degree(3):number(3)
                                    for ii = 1:degree(1)-1
                                        P(1:3,i+ii,j,k) = ((degree(1)-ii)*P(1:3,i,j,k) + ii*P(1:3,i+degree(1),j,k))/degree(1);
                                    end
                                end
                            end
                        end
                        for i = 1:number(1)
                            for j = 1:degree(2):number(2)-degree(2)
                                for k = 1:degree(3):number(3)
                                    for jj = 1:degree(2)-1
                                        P(1:3,i,j+jj,k) = ((degree(2)-jj)*P(1:3,i,j,k) + jj*P(1:3,i,j+degree(2),k))/degree(2);
                                    end
                                end
                            end
                        end
                        for i = 1:number(1)
                            for j = 1:number(2)
                                for k = 1:degree(3):number(3)-degree(3)
                                    for kk = 1:degree(3)-1
                                        P(1:3,i,j,k+kk) = ((degree(3)-kk)*P(1:3,i,j,k) + kk*P(1:3,i,j,k+degree(3)))/degree(3);
                                    end
                                end
                            end
                        end   
                        for i = 1:number(1)
                            for j = 1:number(2)
                                for k = 1:number(3)
                                    P(4,i,j,k) = 1;
                                end
                            end
                        end
                    case 'linear_FEM'
                        P = P(:,1:degree(1):end,1:degree(2):end,1:degree(3):end);
                        P(4,:,:,:) = 1;
                        t{1} = [0, unique(t{1}), 1];
                        t{2} = [0, unique(t{2}), 1];
                        t{3} = [0, unique(t{3}), 1];
                        degree = [1,1,1];
                        number = [length(t{1})-2,...
                                  length(t{2})-2,...
                                  length(t{3})-2];
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
                        number = zeros(1,d_p);
                        t = cell(1,d_p);
                        xi_t = cell(1,d_p);
                        degree = min(nurbs.degree,task.msh.degree);
                        for i = 1:nurbs.d_p
                            t{i} = [0,repelem(unique(nurbs.knots{i}),degree(i)),1];
                            number(i) = numel(t{i})-(degree(i)+1);
                            xi_t{i} = aveknt(t{i}, degree(i)+1);
%                             xi_t{1}{i} = splinesGLL(t{i}, degree(i)+1);
                        end

                        % Evaluate original nurbs at Greville abscissa for interpolation
                        [xi,eta] = ndgrid(xi_t{1},xi_t{2});                        
                        Q = reshape(evaluateNURBSvec(nurbs,[xi(:),eta(:)]).',[nurbs.d,number]);


                        B_xi_arr = zeros(number(1),degree(1)+1);
                        i_xi_arr = zeros(number(1),1);
                        parfor i = 1:number(1)
                            i_xi_p = findKnotSpan(number(1), degree(1), xi_t{1}(i), t{1})        
                            B_xi_arr(i,:) = BsplineBasis(i_xi_p, xi_t{1}(i), degree(1), t{1});
                            i_xi_arr(i) = i_xi_p;
                        end

                        B_eta_arr = zeros(number(2),degree(2)+1);
                        j_eta_arr = zeros(number(2),1);
                        parfor j = 1:number(2)
                            j_eta_p = findKnotSpan(number(2), degree(2), xi_t{2}(j), t{2})        
                            B_eta_arr(j,:) = BsplineBasis(j_eta_p, xi_t{2}(j), degree(2), t{2});
                            j_eta_arr(j) = j_eta_p;
                        end

                        index = zeros(prod(number),2);
                        counter = 1;
                        for j = 1:number(2)
                            for i = 1:number(1)   
                                index(counter,:) = [i, j];
                                counter = counter + 1;
                            end
                        end

                        values = zeros(prod(number),prod(degree+1));
                        n_idx = zeros(prod(number),prod(degree+1));
                        m_idx = zeros(prod(number),prod(degree+1));
                        parfor a = 1:prod(number)
                            i = index(a,1);
                            j = index(a,2);

                            i_xi   = i_xi_arr(i);
                            B_xi = B_xi_arr(i,:);
                            j_eta   = j_eta_arr(j);
                            B_eta = B_eta_arr(j,:);
                            m = i + (j-1)*number(1);

                            temp = zeros(1,prod(degree+1));
                            temp2 = zeros(1,prod(degree+1));
                            counter = 1;
                            for jj = 1:degree(2)+1
                                j_tilde = j_eta - degree(2) + jj - 1;
                                for ii = 1:degree(1)+1
                                    i_tilde = i_xi - degree(1) + ii - 1;
                                    n = i_tilde + (j_tilde-1)*number(1);

                                    temp(counter) = B_xi(ii)*B_eta(jj);
                                    temp2(counter) = n;
                                    counter = counter + 1;
                                end
                            end
                            values(a,:) = temp;
                            n_idx(a,:) = temp2;
                            m_idx(a,:) = ones(1,prod(degree+1))*m
                        end
                        A = sparse(m_idx,n_idx,values);

                        Q = reshape(Q,3,prod(number));
                        P = ones([nurbs.d+1,number]);
                        P_tilde = zeros(size(Q));
                        P_tilde(1,:) = (A\Q(1,:)')';
                        P_tilde(2,:) = (A\Q(2,:)')';
                        P_tilde(3,:) = (A\Q(3,:)')';
                        P(1:3,:,:,:) = reshape(P_tilde,3,number(1),number(2));
                    case 'h_FEM'
                        number(1) = nurbs.number(1);
                        number(2) = nurbs.number(2);

                        degree(1) = nurbs.degree(1);
                        degree(2) = nurbs.degree(2);

                        for i = 1:degree(1):number(1)-degree(1)
                            for j = 1:degree(2):number(2)
                                for ii = 1:degree(1)-1
                                    P(1:3,i+ii,j) = ((degree(1)-ii)*P(1:3,i,j) + ii*P(1:3,i+degree(1),j))/degree(1);
                                end
                            end
                        end
                        for i = 1:number(1)
                            for j = 1:degree(2):number(2)-degree(2)
                                for jj = 1:degree(2)-1
                                    P(1:3,i,j+jj) = ((degree(2)-jj)*P(1:3,i,j) + jj*P(1:3,i,j+degree(2)))/degree(2);
                                end
                            end
                        end 
                        for i = 1:number(1)
                            for j = 1:number(2)
                                P(4,i,j) = 1;
                            end
                        end
                    case 'linear_FEM'
                        P = P(:,1:degree(1):end,1:degree(2):end);
                        P(4,:,:) = 1;
                        t{1} = [0, unique(t{1}), 1];
                        t{2} = [0, unique(t{2}), 1];
                        degree = [1,1];
                        number = [length(t{1})-2,...
                                  length(t{2})-2];

                end
        end
        nurbs.degree = degree;
        nurbs.number = number;
        nurbs.knots = t;
        nurbs.coeffs = P;
        nurbsPatches{patch} = nurbs;
    end
    task.varCol{i_varCol}.nurbs = nurbsPatches;
end
