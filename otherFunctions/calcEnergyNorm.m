function relError = calcEnergyNorm(varColCell, U_cell, options)

formulation = varColCell{1}.formulation;
unconj = strcmp(formulation, 'PGU') || strcmp(formulation,'BGU');

noDomains = length(varColCell);
nodes = cell(noDomains,1);
factors = cell(noDomains,1);
u_hs = cell(noDomains,1);
du_hs = cell(noDomains,1);
for i = 1:noDomains
    U = U_cell{i};
    varCol = varColCell{i};
    %% Extract all needed data from varCol
    Xi = varCol.nurbs.knots{1};
    Eta = varCol.nurbs.knots{2};
    Zeta = varCol.nurbs.knots{3};
    p_xi = varCol.nurbs.degree(1);
    p_eta = varCol.nurbs.degree(2);
    p_zeta = varCol.nurbs.degree(3);

    index = varCol.index;
    noElems = varCol.noElems;
    element = varCol.element;
    elRangeXi = varCol.elRangeXi;
    elRangeEta = varCol.elRangeEta;
    elRangeZeta = varCol.elRangeZeta;
    weights = varCol.weights;
    controlPts = varCol.controlPts;
    
    noCtrlPts = varCol.noCtrlPts;

    %% Preallocation and initiallizations
    if noElems < 600
%         noQuadPts = 30;
        noQuadPts = 16;
    else
        noQuadPts = 10;
    end
    [W3D,Q3D] = gaussianQuadNURBS(noQuadPts,noQuadPts,noQuadPts); 

    if mod(i,2) == 0
        noComponents = 3;
        noComponentsDeriv = 6;
        
        Ux = U(1:3:3*noCtrlPts);
        Uy = U(2:3:3*noCtrlPts);
        Uz = U(3:3:3*noCtrlPts);
    else
        Ux = NaN;
        Uy = NaN;
        Uz = NaN;
        noComponents = 1;
        noComponentsDeriv = 3;
    end
    u_h = zeros(size(W3D,1),noElems,noComponents);
    du_h = zeros(size(W3D,1),noElems,noComponentsDeriv);
    fact = zeros(size(W3D,1),noElems);
    points = zeros(size(W3D,1),noElems,3);

    parfor e = 1:noElems
% keyboard
%     for e = 1:noElems
        idXi = index(e,1);
        idEta = index(e,2);
        idZeta = index(e,3);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);
        Zeta_e = elRangeZeta(idZeta,:);

        J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        u_h_temp = zeros(size(W3D,1),noComponents);
        du_h_temp = zeros(size(W3D,1),noComponentsDeriv);
        fact_temp = zeros(size(W3D,1),1);
        points_temp = zeros(size(W3D,1),3);
        sctr = element(e,:);
        if mod(i,2) == 0
            Uxyz = [Ux(sctr) Uy(sctr) Uz(sctr)];
        else
            U_sctr = U(sctr,:);
        end
        pts = controlPts(sctr,:);

        for gp = 1:size(W3D,1)
            pt = Q3D(gp,:);
            wt = W3D(gp);

            xi   = parent2ParametricSpace(Xi_e,  pt(1));
            eta  = parent2ParametricSpace(Eta_e, pt(2));
            zeta = parent2ParametricSpace(Zeta_e,pt(3));

            [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);

            J = pts'*[dRdxi' dRdeta' dRdzeta'];
            J_1 = det(J);
            X = R*pts;
            dRdX = J'\[dRdxi; dRdeta; dRdzeta];
            if mod(i,2) == 0
                u_h_temp(gp,:) = R*Uxyz;
                
                dUdx = dRdX(1,:)*Uxyz;
                dUdy = dRdX(2,:)*Uxyz;
                dUdz = dRdX(3,:)*Uxyz;
                
                du_h_temp(gp,:) = calculateStrainVector(dUdx, dUdy, dUdz);
            else
                u_h_temp(gp,:) = R*U_sctr;
                du_h_temp(gp,:) = dRdX*U_sctr;
            end
            fact_temp(gp) = J_1 * J_2 * wt;
            points_temp(gp,:) = X;
        end
        u_h(:,e,:) = u_h_temp;
        du_h(:,e,:) = du_h_temp;
        fact(:,e) = fact_temp;
        points(:,e,:) = points_temp;
    end
    u_hs{i} = reshape(u_h, size(u_h,1)*size(u_h,2),noComponents);
    du_hs{i} = reshape(du_h, size(du_h,1)*size(du_h,2),noComponentsDeriv);
    factors{i} = reshape(fact, size(fact,1)*size(fact,2),1);
    nodes{i} = reshape(points, size(points,1)*size(points,2),3);
end

if strcmp(varColCell{1}.applyLoad, 'radialPulsation')
    data.p = varColCell{1}.analytic(nodes{1});
    dp = varColCell{1}.gAnalytic(nodes{1});
    data.dpdx = dp(:,1);
    data.dpdy = dp(:,2);
    data.dpdz = dp(:,3);
else
    data = e3Dss(nodes,options);
end
% keyboard
c_f = options.c_f;
omega = options.omega;
rho_f = options.rho_f;
if noDomains > 1
    rho_s = options.rho_s;
end
k = omega./c_f;
Error = 0;
normalization = 0;
m = 1;
for i = 1:noDomains  
    if mod(i,2) == 0
        sigma = [data(m).sigma_xx,data(m).sigma_yy,data(m).sigma_zz,data(m).sigma_yz,data(m).sigma_xz,data(m).sigma_xy];
        u = [data(m).u_x,data(m).u_y,data(m).u_z];
        C = varColCell{i}.C;
        strain_vec = (C\sigma.').';
        
        strain_error = strain_vec-du_hs{i};
        if unconj
            a_e_h = sum((strain_error*C).*strain_error,2);
            a_u = sum((strain_vec*C).*strain_vec,2);
            u_e = sum(u.*u,2);
            u_e2 = sum((u-u_hs{i}).*(u-u_hs{i}),2);
        else
            a_e_h = sum((strain_error*C).*conj(strain_error),2);
            a_u = sum((strain_vec*C).*conj(strain_vec),2);
            u_e = sum(u.*conj(u),2);
            u_e2 = sum((u-u_hs{i}).*conj(u-u_hs{i}),2);
        end
        
        Error = Error + sum((a_e_h - rho_s(m)*omega^2*u_e2).*factors{i});
        normalization = normalization + sum((a_u - rho_s(m)*omega^2*u_e).*factors{i});
        m = m + 1;
    else
        dp = [data(m).dpdx, data(m).dpdy, data(m).dpdz];
        p = data(m).p;
        if unconj
            Error = Error + 1/(rho_f(m)*omega^2)*sum((sum((dp-du_hs{i}).*(dp-du_hs{i}),2)-k(m)^2*(p-u_hs{i}).^2).*factors{i});
            normalization = normalization + 1/(rho_f(m)*omega^2)*sum((sum(dp.*dp,2)-k(m)^2*p.^2).*factors{i});
        else
            Error = Error + 1/(rho_f(m)*omega^2)*sum((sum((dp-du_hs{i}).*conj(dp-du_hs{i}),2)-k(m)^2*(p-u_hs{i}).*conj(p-u_hs{i})).*factors{i});
            normalization = normalization + 1/(rho_f(m)*omega^2)*sum((sum(dp.*conj(dp),2)-k(m)^2*p.*conj(p)).*factors{i});
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate contribution from interfaces

normalsGlob = cell(noDomains-1,1);
factors = cell(noDomains-1,1);
u_hs = cell(noDomains-1,1);
p_hs = cell(noDomains-1,1);
u = cell(noDomains-1,1);
p = cell(noDomains-1,1);
for i = 1:noDomains-1
    if mod(i,2)
        U = U_cell{i+1};
        P = U_cell{i};
        varCol = varColCell{i+1};
    else
        U = U_cell{i};
        P = U_cell{i+1};
        varCol = varColCell{i};
    end
    
    %% Extract all needed data from varCol
    Xi = varCol.nurbs.knots{1};
    Eta = varCol.nurbs.knots{2};
    p_xi = varCol.nurbs.degree(1);
    p_eta = varCol.nurbs.degree(2);
    n_xi = varCol.nurbs.number(1);
    n_eta = varCol.nurbs.number(2);
    n_zeta = varCol.nurbs.number(3);

    
    elRangeXi = varCol.elRangeXi;
    elRangeEta = varCol.elRangeEta;
    weights = varCol.weights;
    controlPts = varCol.controlPts;
    
    noCtrlPts = varCol.noCtrlPts;
    gluedNodes = varCol.gluedNodes;

    Ux = U(1:3:3*noCtrlPts);
    Uy = U(2:3:3*noCtrlPts);
    Uz = U(3:3:3*noCtrlPts);
    
    [solidXiEtaMesh, solidIndexXiEta, noElems] = generateIGA2DMesh(Xi, Eta, p_xi, p_eta, n_xi, n_eta);
    
    if noElems < 600
        noQuadPts = 16;
    else
        noQuadPts = 10;
    end
    [W2D,Q2D] = gaussianQuadNURBS(noQuadPts,noQuadPts); 
    u_h = zeros(size(W2D,1),noElems,3);
    p_h = zeros(size(W2D,1),noElems);
    fact = zeros(size(W2D,1),noElems);
    points = zeros(size(W2D,1),noElems,3);
    normals = zeros(size(W2D,1),noElems,3);
    

    if mod(i,2)
        solidNodes = zeros(1,n_xi*n_eta);
        counter = 1;
        for jj = 1:n_eta
            for ii = 1:n_xi
                solidNodes(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(jj-1) + ii;
                counter = counter + 1;
            end
        end
        fluidNodes = zeros(1,n_xi*n_eta);
        counter = 1;
        for jj = 1:n_eta
            for ii = 1:n_xi
                fluidNodes(counter) = n_xi*(jj-1) + ii;
                counter = counter + 1;
            end
        end
    else
        solidNodes = zeros(1,n_xi*n_eta);
        counter = 1;
        for jj = 1:n_eta
            for ii = 1:n_xi
                solidNodes(counter) = n_xi*(jj-1) + ii;
                counter = counter + 1;
            end
        end
        varCol_i = varColCell{i+1};
        n_xi = varCol_i.nurbs.number(1);
        n_eta = varCol_i.nurbs.number(2);
        n_zeta = varCol_i.nurbs.number(3);
        fluidNodes = zeros(1,n_xi*n_eta);
        counter = 1;
        for jj = 1:n_eta
            for ii = 1:n_xi
                fluidNodes(counter) = (n_eta*n_xi)*(n_zeta-1) + n_xi*(jj-1) + ii;
                counter = counter + 1;
            end
        end
    end


    % Glue nodes in 2D mesh
    for ii = 1:length(gluedNodes)
        parentIdx = gluedNodes{ii}(1);
        for jj = 2:length(gluedNodes{ii})
            indices = (solidXiEtaMesh == gluedNodes{ii}(jj));
            solidXiEtaMesh(indices) = parentIdx;
        end
    end
    

    parfor e = 1:noElems
% keyboard
%     for e = 1:noElems
        idXi = solidIndexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
        idEta = solidIndexXiEta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        solidSctrXiEta = solidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
        
        fluidSctrXiEta = fluidNodes(solidXiEtaMesh(e,:));          %  element scatter vector
        pts = controlPts(solidSctrXiEta,:);


        u_h_temp = zeros(size(W2D,1),3);
        p_h_temp = zeros(size(W2D,1),1);
        fact_temp = zeros(size(W2D,1),1);
        points_temp = zeros(size(W2D,1),3);
        normals_temp = zeros(size(W2D,1),3);
        
        U_sctr = [Ux(solidSctrXiEta) Uy(solidSctrXiEta) Uz(solidSctrXiEta)];
        P_sctr = P(fluidSctrXiEta);

        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta = parent2ParametricSpace(Eta_e,pt(2));

            [R, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights(solidNodes));

            J = pts'*[dRxi' dRdeta'];
            crossProd = cross(J(:,1), J(:,2));  % pointing outwards if sphere
            J_1 = norm(crossProd);
            n = crossProd/J_1;

            X = R*pts;
            u_h_temp(gp,:) = R*U_sctr;
            p_h_temp(gp) = R*P_sctr;
            fact_temp(gp) = J_1 * J_2 * wt;
            points_temp(gp,:) = X;
            normals_temp(gp,:) = n;
        end
        u_h(:,e,:) = u_h_temp;
        p_h(:,e,:) = p_h_temp;
        fact(:,e) = fact_temp;
        points(:,e,:) = points_temp;
        normals(:,e,:) = normals_temp;
    end
    u_hs{i} = reshape(u_h, size(u_h,1)*size(u_h,2),3);
    p_hs{i} = reshape(p_h, size(p_h,1)*size(p_h,2),1);
    factors{i} = reshape(fact, size(fact,1)*size(fact,2),1);
    
    nodes = cell(noDomains,1);
    nodes{i} = reshape(points, size(points,1)*size(points,2),3);
    nodes{i+1} = nodes{i};
    data = e3Dss(nodes,options);
    if mod(i,2)
        u{i} = [data(i).u_x,data(i).u_y,data(i).u_z];
        p{i} = data(i).p;
    else
        u{i} = [data(i-1).u_x,data(i-1).u_y,data(i-1).u_z];
        p{i} = data(i).p;
    end
    
    normalsGlob{i} = reshape(normals, size(normals,1)*size(normals,2),3);
end

m = 1;
for i = 1:noDomains-1
    u_r_error = sum((u{i}-u_hs{i}).*normalsGlob{i},2);
    p_error = p{i}-p_hs{i};
    
    u_r = sum(u{i}.*normalsGlob{i},2);
    if unconj % terms not similar otherwise
        Error = Error + 2*sum(u_r_error.*p_error.*factors{i});
        normalization = normalization + 2*sum(u_r.*p{i}.*factors{i});
    else
        Error = Error + sum((u_r_error.*conj(p_error) + p_error.*conj(u_r_error)).*factors{i});
        normalization = normalization + sum((u_r.*conj(p{i})+p{i}.*conj(u_r)).*factors{i});        
    end
    m = m + 1;
end

relError = 100*sqrt(abs(Error)/abs(normalization));
