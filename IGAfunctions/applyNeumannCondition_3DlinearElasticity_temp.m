if ~exist('F','var')
    F = zeros(noDofs,1);        % external force vector
elseif length(F) ~= noDofs
    F = zeros(noDofs,1);        % external force vector
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
    for e = 1:noElemsXiEta
%     parfor e = 1:noElemsXiEta
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
%             switch coordinateSystem
%                 case 0
%                     StressMatrix6 = stressMatrix(stress_u(x,y,z));
%                 case 1                    
%                     radius = sqrt(x^2+y^2);
%                     theta = atan2(y,x);
%                     StressMatrix6 = stressMatrix(stress_u(radius,theta,z));
%             end
            normal = crossProd/norm(crossProd);
            keyboard
%             traction = StressMatrix6*normal;

%         F1_e = F1_e + P_inc(x,y,z)*kron(R_fun',normal)*norm(crossProd) * J_2 * wt;
            f_e = f_e + P_inc(x,y,z)*kron(normal,R_fun') * norm(crossProd) * J_2 * wt;  
%             f_e = f_e + [traction(1)*R_fun'; traction(2)*R_fun'; traction(3)*R_fun'] * norm(crossProd) * J_2 * wt;  
        end   
        indices(:,e) = sctrXiEtaVec';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end
