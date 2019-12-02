if ~exist('F','var')
    F = zeros(noDofs,1);        % external force vector
elseif length(F) ~= noDofs
    F = zeros(noDofs,1);        % external force vector
end
n = varCol{2}.nurbs.number(1);
m = varCol{2}.nurbs.number(2);
l = varCol{2}.nurbs.number(3);
weights = varCol{2}.weights;
controlPts = varCol{2}.controlPts;
p = varCol{2}.nurbs.degree(1);
q = varCol{2}.nurbs.degree(2);
r = varCol{2}.nurbs.degree(3);
Xi = varCol{2}.nurbs.knots{1};
Eta = varCol{2}.nurbs.knots{2};
Zeta = varCol{2}.nurbs.knots{3};
gluedNodes = varCol{2}.gluedNodes;
elRangeXi = varCol{2}.elRangeXi;
elRangeEta = varCol{2}.elRangeEta;
elRangeZeta = varCol{2}.elRangeZeta;
noCtrlPts = varCol{2}.noCtrlPts;
noElems = varCol{2}.noElems;
index = varCol{2}.index;
element = varCol{2}.element;
nurbs = varCol{2}.nurbs;

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
    Fvalues = zeros(3*n_en,noElemsEtaZeta);
    indices = zeros(3*n_en,noElemsEtaZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(q+1,r+1); 
    parfor e = 1:noElemsEtaZeta
%     for e = 1:noElemsEtaZeta
        idEta = indexEtaZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexEtaZeta(e,2);

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrEtaZeta = xi0Nodes(EtaZetaMesh(e,:));          %  element scatter vector
        sctrEtaZetaVec = [sctrEtaZeta sctrEtaZeta+noCtrlPts sctrEtaZeta+2*noCtrlPts];
        n_en = length(sctrEtaZeta);
        pts = controlPts(sctrEtaZeta,:);
        f_e = zeros(3*n_en,1);
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
                    StressMatrix1 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix1 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix1 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = -crossProd/norm(crossProd);
            traction = StressMatrix1*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrEtaZetaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
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
    Fvalues = zeros(3*n_en,noElemsEtaZeta);
    indices = zeros(3*n_en,noElemsEtaZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(q+1,r+1); 
    parfor e = 1:noElemsEtaZeta
%     for e = 1:noElemsEtaZeta
        idEta = indexEtaZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexEtaZeta(e,2);

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrEtaZeta = xi1Nodes(EtaZetaMesh(e,:));          %  element scatter vector
        sctrEtaZetaVec = [sctrEtaZeta sctrEtaZeta+noCtrlPts sctrEtaZeta+2*noCtrlPts];
        n_en = length(sctrEtaZeta);
        pts = controlPts(sctrEtaZeta,:);
        f_e = zeros(3*n_en,1);
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
            StressMatrix2 = zeros(6);
            switch coordinateSystem
                case 0
                    StressMatrix2 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix2 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix2 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = crossProd/norm(crossProd);
            traction = StressMatrix2*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrEtaZetaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
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
    Fvalues = zeros(3*n_en,noElemsXiZeta);
    indices = zeros(3*n_en,noElemsXiZeta);
    
    [W2D,Q2D] = gaussianQuadNURBS(p+1,r+1); 
    parfor e = 1:noElemsXiZeta
%     for e = 1:noElemsXiZeta
        idXi = indexXiZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexXiZeta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrXiZeta = eta0Nodes(XiZetaMesh(e,:));          %  element scatter vector
        sctrXiZetaVec = [sctrXiZeta sctrXiZeta+noCtrlPts sctrXiZeta+2*noCtrlPts];

        n_en = length(sctrXiZeta);
        pts = controlPts(sctrXiZeta,:);
        f_e = zeros(3*n_en,1);
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
                    StressMatrix3 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix3 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix3 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = crossProd/norm(crossProd);
            traction = StressMatrix3*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXiZetaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
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
    Fvalues = zeros(3*n_en,noElemsXiZeta);
    indices = zeros(3*n_en,noElemsXiZeta);
 
    [W2D,Q2D] = gaussianQuadNURBS(p+1,r+1); 
    parfor e = 1:noElemsXiZeta
%     for e = 1:noElemsXiZeta
        idXi = indexXiZeta(e,1);   % the index matrix is made in generateIGA3DMesh
        idZeta = indexXiZeta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctrXiZeta = eta1Nodes(XiZetaMesh(e,:));          %  element scatter vector
        sctrXiZetaVec = [sctrXiZeta sctrXiZeta+noCtrlPts sctrXiZeta+2*noCtrlPts];

        n_en = length(sctrXiZeta);
        pts = controlPts(sctrXiZeta,:);
        f_e = zeros(3*n_en,1);
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
                    StressMatrix4 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix4 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix4 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = -crossProd/norm(crossProd);
            traction = StressMatrix4*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXiZetaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
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
    Fvalues = zeros(3*n_en,noElemsXiEta);
    indices = zeros(3*n_en,noElemsXiEta);
    
    [W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
%     parfor e = 1:noElemsXiEta
    for e = 1:noElemsXiEta
        idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
        idEta = indexXiEta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector
        sctrXiEtaVec = [sctrXiEta sctrXiEta+noCtrlPts sctrXiEta+2*noCtrlPts];

        n_en = length(sctrXiEta);
        pts = controlPts(sctrXiEta,:);
        f_e = zeros(3*n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta = parent2ParametricSpace(Eta_e,pt(2));
            
            [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

            J = pts'*[dRdxi' dRdeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    StressMatrix5 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix5 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix5 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = -crossProd/norm(crossProd);
            traction = StressMatrix5*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXiEtaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
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
    Fvalues = zeros(3*n_en,noElemsXiEta);
    indices = zeros(3*n_en,noElemsXiEta);
    
    [W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
    
%     parfor e = 1:noElemsXiEta
    for e = 1:noElemsXiEta
        idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
        idEta = indexXiEta(e,2);

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
        Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector
        sctrXiEtaVec = [sctrXiEta sctrXiEta+noCtrlPts sctrXiEta+2*noCtrlPts];

        n_en = length(sctrXiEta);
        pts = controlPts(sctrXiEta,:);
        f_e = zeros(3*n_en,1);
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta = parent2ParametricSpace(Eta_e,pt(2));
            
            [R_fun, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

            J = pts'*[dRdxi' dRdeta'];
            crossProd = cross(J(:,1), J(:,2)); 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            z = v(3);
            switch coordinateSystem
                case 0
                    StressMatrix6 = stressMatrix(stress_u(x,y,z));
                case 1                    
                    radius = sqrt(x^2+y^2);
                    theta = atan2(y,x);
                    StressMatrix6 = stressMatrix(stress_u(radius,theta,z));
                case 2                    
                    phi = atan2(y,x);
                    radius = sqrt(x^2+y^2+z^2);
                    theta = acos(z/radius);
                    StressMatrix6 = stressMatrix(stress_u(radius,theta,phi));
            end
            normal = crossProd/norm(crossProd);
            traction = StressMatrix6*normal;

            f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end   
        indices(:,e) = sctrXiEtaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end
