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

                    case {'hp_FEM','sub_IGA'}
                        number = zeros(1,d_p);
                        t = cell(1,d_p);
                        xi_t = cell(1,d_p);
                        degree = min(nurbs.degree,task.msh.degree);
                        for i = 1:nurbs.d_p
                            switch task.misc.coreMethod
                                case 'hp_FEM'
                                    t{i} = [0,repelem(unique(nurbs.knots{i}),degree(i)),1];
                                case 'sub_IGA'
                                    % Lower order while maintaining continuity
                                    t{i} = nurbs.knots{i};
                                    uniqueKnots = unique(t{i});
                                    for j = 1:(nurbs.degree(i) - task.msh.degree)
                                        t{i} = setdiffUnique(t{i},uniqueKnots);
                                        t{i} = setdiffUnique(t{i},uniqueKnots);
                                        t{i} = sort([t{i},uniqueKnots]);
                                    end
                            end
                            number(i) = numel(t{i})-(degree(i)+1);
                            xi_t{i} = aveknt(t{i}, degree(i)+1);
%                             xi_t{1}{i} = splinesGLL(t{i}, degree(i)+1);
                        end

                        % Evaluate original nurbs at Greville abscissa for interpolation
                        [xi,eta,zeta] = ndgrid(xi_t{1},xi_t{2},xi_t{3});                        
                        Q = reshape(evaluateNURBS(nurbs,[xi(:),eta(:),zeta(:)]).',[nurbs.d,number]);

                        i_arr = cell(1,d_p);
                        B_arr = cell(1,d_p);
                        for i = 1:d_p
                            i_arr{i} = findKnotSpanVec(number(i), degree(i), xi_t{i}.', t{i});   
                            B_arr{i} = BsplineBasisVec(i_arr{i}, xi_t{i}.', degree(i), t{i});
                        end

                        [arr1,arr2,arr3] = ndgrid(1:number(1),1:number(2),1:number(3));
                        index = [arr1(:),arr2(:),arr3(:)];

                        values = zeros(prod(number),prod(degree+1));
                        n_idx = zeros(prod(number),prod(degree+1));
                        m_idx = zeros(prod(number),prod(degree+1));
                        ii = 1:degree(1)+1;
                        jj = 1:degree(2)+1;
                        ll = 1:degree(3)+1;
                        parfor a = 1:prod(number)
                            i = index(a,1);
                            j = index(a,2);
                            l = index(a,3);

                            i_tilde = i_arr{1}(i) - degree(1) + ii - 1;
                            j_tilde = i_arr{2}(j) - degree(2) + jj - 1;
                            l_tilde = i_arr{3}(l) - degree(3) + ll - 1;
                            values(a,:) = reshape(B_arr{1}(i,:).'.*B_arr{2}(j,:).*reshape(B_arr{3}(l,:),1,1,[]),1,[]);
                            n_idx(a,:) = reshape(i_tilde.' + (j_tilde-1)*number(1) + (reshape(l_tilde,1,1,[])-1)*number(1)*number(2),1,[]);
                            m_idx(a,:) = ones(1,prod(degree+1))*(i + (j-1)*number(1) + (l-1)*number(1)*number(2));
                        end
                        A = sparse(m_idx,n_idx,values);

                        Q = reshape(Q,3,prod(number));
                        P_tilde = (A\Q.').';
                        P = ones([nurbs.d+1,number]);
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

                    case {'hp_FEM','sub_IGA'}
                        number = zeros(1,d_p);
                        t = cell(1,d_p);
                        xi_t = cell(1,d_p);
                        degree = min(nurbs.degree,task.msh.degree);
                        for i = 1:nurbs.d_p
                            switch task.misc.coreMethod
                                case 'hp_FEM'
                                    t{i} = [0,repelem(unique(nurbs.knots{i}),degree(i)),1];
                                case 'sub_IGA'
                                    % Lower order while maintaining continuity
                                    uniqueKnots = unique(nurbs.knots{i});
                                    t{i} = setdiffUnique(nurbs.knots{i},uniqueKnots);
                                    t{i} = setdiffUnique(t{i},uniqueKnots);
                                    t{i} = sort([t{i},uniqueKnots]);
                            end
                            number(i) = numel(t{i})-(degree(i)+1);
                            xi_t{i} = aveknt(t{i}, degree(i)+1);
%                             xi_t{1}{i} = splinesGLL(t{i}, degree(i)+1);
                        end

                        % Evaluate original nurbs at Greville abscissa for interpolation
                        [xi,eta] = ndgrid(xi_t{1},xi_t{2});                        
                        Q = reshape(evaluateNURBS(nurbs,[xi(:),eta(:)]).',[nurbs.d,number]);

                        i_arr = cell(1,d_p);
                        B_arr = cell(1,d_p);
                        for i = 1:d_p
                            i_arr{i} = findKnotSpanVec(number(i), degree(i), xi_t{i}.', t{i});   
                            B_arr{i} = BsplineBasisVec(i_arr{i}, xi_t{i}.', degree(i), t{i});
                        end

                        [arr1,arr2] = ndgrid(1:number(1),1:number(2));
                        index = [arr1(:),arr2(:)];


                        values = zeros(prod(number),prod(degree+1));
                        n_idx = zeros(prod(number),prod(degree+1));
                        m_idx = zeros(prod(number),prod(degree+1));
                        ii = 1:degree(1)+1;
                        jj = 1:degree(2)+1;
                        parfor a = 1:prod(number)
                            i = index(a,1);
                            j = index(a,2);

                            i_tilde = i_arr{1}(i) - degree(1) + ii - 1;
                            j_tilde = i_arr{2}(j) - degree(2) + jj - 1;
                            values(a,:) = reshape(B_arr{1}(i,:).'.*B_arr{2}(j,:),1,[]);
                            n_idx(a,:) = reshape(i_tilde.' + (j_tilde-1)*number(1),1,[]);
                            m_idx(a,:) = ones(1,prod(degree+1))*(i + (j-1)*number(1));
                        end
                        A = sparse(m_idx,n_idx,values);

                        Q = reshape(Q,3,prod(number));
                        P_tilde = (A\Q.').';
                        P = ones([nurbs.d+1,number]);
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
