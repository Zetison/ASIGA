function F = applyHWBC_2D(varCol,alpha_s_arr,x_0,k_wn,P_0)

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

elRangeXi = varCol.elRangeXi;

Xi = varCol.nurbs.knots{1};
p = varCol.nurbs.degree(1);
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);

gluedNodes = varCol.gluedNodes;

noDofs_tot = varCol.noDofs_tot;

weights = varCol.weights;
controlPts = varCol.controlPts;

rightHandedOrientation = findOrientation_2D(varCol);


% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs_tot,length(alpha_s_arr));        % external force vector


% Apply Neumann boundary condition at surface where eta = 0
eta0Nodes = zeros(1,n);
counter = 1;
for i = 1:n
    eta0Nodes(counter) = i;
    counter = counter + 1;
end

[XiMesh, indexXi, noElemsXi] = generateIGA1DMesh(Xi, p, n);

% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (XiMesh == gluedNodes{i}(j));
        XiMesh(indices) = parentIdx;
    end
end

n_en = p+1;
Fvalues = zeros(n_en,noElemsXi,length(alpha_s_arr));
indices = zeros(n_en,noElemsXi);

[W1D,Q1D] = gaussianQuadNURBS(p+1); 

% parfor e = 1:noElemsXi

for e = 1:noElemsXi
    idXi = indexXi(e,1);   % the index matrix is made in generateIGA3DMesh

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]

    J_2 = 0.5*(Xi_e(2)-Xi_e(1));

    sctrXiEta = eta0Nodes(XiMesh(e,:));          %  element scatter vector

    n_en = length(sctrXiEta);
    pts = controlPts(sctrXiEta,:);
    f_e = zeros(n_en,length(alpha_s_arr));
    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));

        [R_fun, dRxi] = NURBS1DBasis(xi, p, Xi, weights(eta0Nodes));

        J = pts'*dRxi';
        crossProd = J'; 
        

        v = R_fun*pts;
        
        tangent = crossProd/norm(crossProd);
        normal = [-tangent(2), tangent(1)];
        k_vec = -[k_wn*cos(alpha_s_arr);
                  k_wn*sin(alpha_s_arr)];
              
        P_inc = P_0*exp(-1i*dot2(k_vec,x_0)).*exp(1i*dot2(k_vec,v'));
        dP_inc = 1i*dot2(k_vec,normal).*P_inc;
        deriv = -dP_inc;
        
        f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctrXiEta';
    Fvalues(:,e,:) = f_e;
end

for alpha_s_Nr = 1:length(alpha_s_arr)
    F(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
end

% 
% 
% % Apply Neumann boundary condition at surface where eta = 1
% eta1Nodes = zeros(1,n);
% counter = 1;
% for i = 1:n
%     eta1Nodes(counter) = n*(m-1) + i;
%     counter = counter + 1;
% end
% 
% [XiMesh, indexXi, noElemsXi] = generateIGA1DMesh(Xi, p, n);
% 
% % Glue nodes in 2D mesh
% for i = 1:length(gluedNodes)
%     parentIdx = gluedNodes{i}(1);
%     for j = 2:length(gluedNodes{i})
%         indices = (XiMesh == gluedNodes{i}(j));
%         XiMesh(indices) = parentIdx;
%     end
% end
% 
% n_en = p+1;
% Fvalues = zeros(n_en,noElemsXi,length(alpha_s_arr));
% indices = zeros(n_en,noElemsXi);
% 
% [W1D,Q1D] = gaussianQuadNURBS(p+1); 
% % parfor e = 1:noElemsXi
% for e = 1:noElemsXi
%     idXi = indexXi(e,1);   % the index matrix is made in generateIGA3DMesh
% 
%     Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
% 
%     J_2 = 0.5*(Xi_e(2)-Xi_e(1));
% 
%     sctrXiEta = eta1Nodes(XiMesh(e,:));          %  element scatter vector
% 
%     n_en = length(sctrXiEta);
%     pts = controlPts(sctrXiEta,:);
%     f_e = zeros(n_en,length(alpha_s_arr));
%     for gp = 1:size(W1D,1)
%         pt = Q1D(gp,:);
%         wt = W1D(gp);
% 
%         xi  = parent2ParametricSpace(Xi_e, pt(1));
% 
%         [R_fun, dRxi] = NURBS1DBasis(xi, p, Xi, weights(eta1Nodes));
% 
%         J = pts'*dRxi';
%         crossProd = J'; 
%         
% 
%         v = R_fun*pts;
%         x = v(1);
%         y = v(2);
%         
%         tangent = -crossProd/norm(crossProd);
%         normal = [-tangent(2), tangent(1)];
%         k_vec = -[k_wn*cos(alpha_s_arr);
%                   k_wn*sin(alpha_s_arr)];
%               
%         deriv = scatteredPressureOnRigidCylinder2_deriv(x, y, P_0, k_wn, R_o, 1e-13,alpha_s_arr);
%         
% %         deriv = Grad(v).'*normal';
%         
%         f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
%     end    
%     indices(:,e) = sctrXiEta';
%     Fvalues(:,e,:) = f_e;
% end
% 
% for alpha_s_Nr = 1:length(alpha_s_arr)
%     F(:,alpha_s_Nr) = F + vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
% end
