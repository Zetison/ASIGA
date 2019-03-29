% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs,1);        % external force vector


if applyNeumannAtXi0
    % Apply Neumann boundary condition at surface where xi = 0
    xi0Nodes = zeros(1,m*l);
    counter = 1;
    for k = 1:l
        for j = 1:m
            xi0Nodes(counter) = (m*n)*(k-1) + n*(j-1) + 1;
            counter = counter + 1;
        end
    end
    
    [EtaZetaMesh, indexEtaZeta, noElemsEtaZeta, ...
     ~, ~] = generateIGA2DMesh(Eta, Zeta, q, r, m, l);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (EtaZetaMesh == gluedNodes{i}(j));
            EtaZetaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (q+1)*(r+1);
    Fvalues = zeros(n_en,noElemsEtaZeta);
    indices = zeros(n_en,noElemsEtaZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(q+1,r+1); 
    parfor e = 1:noElemsEtaZeta
        idEta = indexEtaZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexEtaZeta(e,2);

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrEtaZeta = xi0Nodes(EtaZetaMesh(e,:));          %  element scatter vector
        n_en = length(sctrEtaZeta);
        pts = controlPts(sctrEtaZeta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            eta  = parent2ParametricSpace(Eta_e, pt(1));
            zeta = parent2ParametricSpace(Zeta_e,pt(2));
            
            [R_fun, dRdeta, dRdzeta] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi0Nodes));

            J = pts'*[dRdeta' dRdzeta'];
            crossProd = cross(J(:,1), J(:,2)); 

            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv1 = g_xi0(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv1 = g_xi0(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv1 = g_xi0(radius, phi, theta);
                case 3
                    normal = -crossProd/norm(crossProd);
                    deriv1 = g_xi0(x,y,z, normal);
            end

            f_e = f_e + R_fun'*deriv1*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrEtaZeta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,length(F));
end

if applyNeumannAtXi1
    % Apply Neumann boundary condition at surface where xi = 1
    xi1Nodes = zeros(1,m*l);
    counter = 1;
    for k = 1:l
        for j = 1:m
            xi1Nodes(counter) = (m*n)*(k-1) + n*(j-1) + n;
            counter = counter + 1;
        end
    end
    
    [EtaZetaMesh, indexEtaZeta, noElemsEtaZeta, ...
     ~, ~] = generateIGA2DMesh(Eta, Zeta, q, r, m, l);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (EtaZetaMesh == gluedNodes{i}(j));
            EtaZetaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (q+1)*(r+1);
    Fvalues = zeros(n_en,noElemsEtaZeta);
    indices = zeros(n_en,noElemsEtaZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(q+1,r+1); 
    parfor e = 1:noElemsEtaZeta
        idEta = indexEtaZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexEtaZeta(e,2);

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrEtaZeta = xi1Nodes(EtaZetaMesh(e,:));          %  element scatter vector
        n_en = length(sctrEtaZeta);
        pts = controlPts(sctrEtaZeta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            eta  = parent2ParametricSpace(Eta_e, pt(1));
            zeta = parent2ParametricSpace(Zeta_e,pt(2));
            
            [R_fun, dRdeta, dRdzeta] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi1Nodes));

            J = pts'*[dRdeta' dRdzeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv2 = g_xi1(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv2 = g_xi1(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv2 = g_xi1(radius, phi, theta);
                case 3
                    normal = crossProd/norm(crossProd);
                    deriv2 = g_xi1(x,y,z, normal);
            end

            f_e = f_e + R_fun'*deriv2*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrEtaZeta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,length(F));
end



if applyNeumannAtEta0
    % Apply Neumann boundary condition at surface where eta = 0
    eta0Nodes = zeros(1,n*l);
    counter = 1;
    for k = 1:l
        for i = 1:n
            eta0Nodes(counter) = (m*n)*(k-1) + i;
            counter = counter + 1;
        end
    end

    [XiZetaMesh, indexXiZeta, noElemsXiZeta, ...
     ~, ~] = generateIGA2DMesh(Xi, Zeta, p, r, n, l);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (XiZetaMesh == gluedNodes{i}(j));
            XiZetaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (p+1)*(r+1);
    Fvalues = zeros(n_en,noElemsXiZeta);
    indices = zeros(n_en,noElemsXiZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(p+1,r+1); 
    parfor e = 1:noElemsXiZeta
        idXi = indexXiZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexXiZeta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrXiZeta = eta0Nodes(XiZetaMesh(e,:));          %  element scatter vector

        n_en = length(sctrXiZeta);
        pts = controlPts(sctrXiZeta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            zeta = parent2ParametricSpace(Zeta_e,pt(2));
            
            [R_fun, dRdxi, dRdzeta] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta0Nodes));

            J = pts'*[dRdxi' dRdzeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv3 = g_eta0(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv3 = g_eta0(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv3 = g_eta0(radius, phi, theta);
                case 3
                    normal = crossProd/norm(crossProd);
                    deriv3 = g_eta0(x,y,z, normal);
            end

            f_e = f_e + R_fun'*deriv3*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXiZeta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,length(F));
end


if applyNeumannAtEta1
    % Apply Neumann boundary condition at surface where eta = 1
    eta1Nodes = zeros(1,n*l);
    counter = 1;
    for k = 1:l
        for i = 1:n
            eta1Nodes(counter) = (m*n)*(k-1) + n*(m-1) + i;
            counter = counter + 1;
        end
    end
    
    [XiZetaMesh, indexXiZeta, noElemsXiZeta, ...
     ~, ~] = generateIGA2DMesh(Xi, Zeta, p, r, n, l);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (XiZetaMesh == gluedNodes{i}(j));
            XiZetaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (p+1)*(r+1);
    Fvalues = zeros(n_en,noElemsXiZeta);
    indices = zeros(n_en,noElemsXiZeta);
 
    [W2D,Q2D] = gaussianQuadNURBS(p+1,r+1); 
    parfor e = 1:noElemsXiZeta
        idXi = indexXiZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexXiZeta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrXiZeta = eta1Nodes(XiZetaMesh(e,:));          %  element scatter vector

        n_en = length(sctrXiZeta);
        pts = controlPts(sctrXiZeta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            zeta = parent2ParametricSpace(Zeta_e,pt(2));
            
            [R_fun, dRdxi, dRdzeta] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta1Nodes));

            J = pts'*[dRdxi' dRdzeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv4 = g_eta1(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv4 = g_eta1(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv4 = g_eta1(radius, phi, theta);
                case 3
                    normal = -crossProd/norm(crossProd);
                    deriv4 = g_eta1(x,y,z, normal);
            end

            f_e = f_e + R_fun'*deriv4*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXiZeta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,length(F));
end



if applyNeumannAtZeta0
    % Apply Neumann boundary condition at surface where zeta = 0
    zeta0Nodes = zeros(1,n*m);
    counter = 1;
    for j = 1:m
        for i = 1:n
            zeta0Nodes(counter) = n*(j-1) + i;
            counter = counter + 1;
        end
    end
    
    [XiEtaMesh, indexXiEta, noElemsXiEta, ...
     ~, ~] = generateIGA2DMesh(Xi, Eta, p, q, n, m);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (XiEtaMesh == gluedNodes{i}(j));
            XiEtaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (p+1)*(q+1);
    Fvalues = zeros(n_en,noElemsXiEta);
    indices = zeros(n_en,noElemsXiEta);

    [W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
    parfor e = 1:noElemsXiEta
%     for e = 1:noElemsXiEta
        idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
        idEta = indexXiEta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector

        n_en = length(sctrXiEta);
        pts = controlPts(sctrXiEta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta = parent2ParametricSpace(Eta_e,pt(2));
            
            [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

            J = pts'*[dRxi' dRdeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv5 = g_zeta0(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv5 = g_zeta0(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv5 = g_zeta0(radius, phi, theta);
                case 3
                    if dot(v,crossProd) > 0
                        normal = -crossProd/norm(crossProd);
                    else
                        normal = crossProd/norm(crossProd);
                    end
                    deriv5 = g_zeta0(x,y,z, normal);                    
            end
            f_e = f_e + R_fun'*deriv5*norm(crossProd) * J_2 * wt;  

        end    
        indices(:,e) = sctrXiEta';
        Fvalues(:,e) = f_e;
    end
    
    F = F + vectorAssembly(Fvalues,indices,length(F));
end

if applyNeumannAtZeta1
    % Apply Neumann boundary condition at surface where zeta = 1
    zeta1Nodes = zeros(1,n*m);
    counter = 1;
    for j = 1:m
        for i = 1:n
            zeta1Nodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
            counter = counter + 1;
        end
    end
    
    [XiEtaMesh, indexXiEta, noElemsXiEta, ...
     ~, ~] = generateIGA2DMesh(Xi, Eta, p, q, n, m);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (XiEtaMesh == gluedNodes{i}(j));
            XiEtaMesh(indices) = parentIdx;
        end
    end
    
    n_en = (p+1)*(q+1);
    Fvalues = zeros(n_en,noElemsXiEta);
    indices = zeros(n_en,noElemsXiEta);
 
    [W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
    parfor e = 1:noElemsXiEta
%     for e = 1:noElemsXiEta
        idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
        idEta = indexXiEta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector

        n_en = length(sctrXiEta);
        pts = controlPts(sctrXiEta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta = parent2ParametricSpace(Eta_e,pt(2));
            
            [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

            J = pts'*[dRxi' dRdeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    deriv6 = g_zeta1(x,y,z);
                case 1
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);
                    deriv6 = g_zeta1(radius,theta,z);
                case 2
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    deriv6 = g_zeta1(radius, phi, theta);
                case 3
                    if dot(v,crossProd) > 0
                        normal = crossProd/norm(crossProd);
                    else
                        normal = -crossProd/norm(crossProd);
                    end
                    deriv6 = g_zeta1(x,y,z, normal);
            end

            f_e = f_e + R_fun'*deriv6*norm(crossProd) * J_2 * wt;  
        end   
        indices(:,e) = sctrXiEta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,length(F));
end
