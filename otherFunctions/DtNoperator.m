function k_e = DtNoperator(W2D,Q2D, R_fun, J_2, pts,p,q,r,Xi,Eta,Zeta,weights)

% Apply DtN operator at surface where zeta = 1
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
Kvalues = zeros(n*m*n_en,noElemsXiEta);
indices = zeros(n*m*n_en,noElemsXiEta);

[W2D,Q2D] = gaussianQuadNURBS(p+1,q+1); 
[W2D2,Q2D2] = gaussianQuadNURBS(p+1,q+1); 
parfor e = 1:noElemsXiEta
    idXi = indexXiEta(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = indexXiEta(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector

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

        f_e = f_e + R_fun'*deriv6*norm(crossProd) * J_2 * wt;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for e2 = 1:noElemsXiEta
            idXi2 = indexXiEta(e2,1);   % the index matrix is made in generateIGA3DMesh
            idEta2 = indexXiEta(e2,2);

            Xi_e2 = elRangeXi(idXi2,:); % [eta_j,eta_j+1]
            Eta_e2 = elRangeEta(idEta2,:); % [zeta_k,zeta_k+1]

            J_22 = 0.25*(Xi_e2(2)-Xi_e2(1))*(Eta_e2(2)-Eta_e2(1));

            sctrXiEta2 = zeta1Nodes(XiEtaMesh(e2,:));          %  element scatter vector

            pts = controlPts(sctrXiEta2,:);
            f_e = zeros(n_en,1);
            for gp2 = 1:size(W2D,1)
                pt2 = Q2D2(gp,:);
                wt2 = W2D2(gp);
                xi2  = parent2ParametricSpace(Xi_e2, pt2(1));
                eta2 = parent2ParametricSpace(Eta_e2,pt2(2));

                [R_fun2, dRxi2, dRdeta2] = NURBS2DBasis(xi2, eta2, p, q, Xi, Eta, weights(zeta1Nodes));
                J2 = pts2'*[dRxi2' dRdeta2'];
                crossProd2 = cross(J2(:,1), J2(:,2)); 

                u_mn_h = zeros(size(R_fun));
                for n_y = 1:5
                    for m_y = -n_y:n_y
                        alpha_mn = sqrt(4*pi/(2*n_y+1)*factorial(n_y+abs(m_y))/(factorial(n_y-abs(m_y))));

                        P = legendre(n,cos_theta);
                        Y_mn = P(abs(m)+1)*exp(1i*m*phi)/alpha_mn;

                        u_mn_h = u_mn_h + R_fun2*conj(Y_mn)*norm(crossProd2)*J_22*wt2;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end   
    indices(:,e) = sctrXiEta';
    Kvalues(:,e) = f_e;
end
F = F + vectorAssembly(Kvalues,indices,noDofs);

    