varCol = varCol_fluid_o;
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
gluedNodes = varCol.gluedNodes;

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;

noCtrlPts = varCol.noCtrlPts;
noDofs = varCol.noDofs;

if ~exist('F','var')
    F = zeros(noDofs,1);        % external force vector
elseif length(F) ~= noDofs
    F = zeros(noDofs,1);        % external force vector
end


if applyNeumannAtXi0
    % Apply Neumann boundary condition at surface where xi = 0
    xi0Nodes = zeros(1,m);
    counter = 1;
    for j = 1:m
        xi0Nodes(counter) = n*(j-1) + 1;
        counter = counter + 1;
    end
    
    [EtaMesh, indexEta, noElemsEta] = generateIGA1DMesh(Eta, q, m);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (EtaMesh == gluedNodes{i}(j));
            EtaMesh(indices) = parentIdx;
        end
    end
    
    n_en = q+1;
    Fvalues = zeros(n_en,noElemsEta);
    indices = zeros(n_en,noElemsEta);
    
    [W1D,Q1D] = gaussianQuadNURBS(q+1); 
    for e = 1:noElemsEta
        idEta = indexEta(e,1);   % the index matrix is made in generateIGA3DMesh

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]

        J_2 = 0.5*(Eta_e(2)-Eta_e(1));

        sctrEta = xi0Nodes(EtaMesh(e,:));          %  element scatter vector
        n_en = length(sctrEta);
        pts = controlPts(sctrEta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W1D,1)
            pt = Q1D(gp,:);
            wt = W1D(gp);

            eta  = parent2ParametricSpace(Eta_e, pt(1));
            
            [R_fun, dRdeta] = NURBS1DBasis(eta, q, Eta, weights(xi0Nodes));

            J = pts'*dRdeta';
            crossProd = J'; 

            v = R_fun*pts;
            tangent = crossProd/norm(crossProd);
            normal = [-tangent(2), tangent(1)];
            deriv = Grad(v).'*normal';

            f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
            
        end    
        indices(:,e) = sctrEta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end

if applyNeumannAtXi1
    % Apply Neumann boundary condition at surface where xi = 1
    xi1Nodes = zeros(1,m);
    counter = 1;
    for j = 1:m
        xi1Nodes(counter) = n*(j-1) + n;
        counter = counter + 1;
    end
    
    [EtaMesh, indexEta, noElemsEta] = generateIGA1DMesh(Eta, q, m);
 
    % Glue nodes in 2D mesh
    for i = 1:length(gluedNodes)
        parentIdx = gluedNodes{i}(1);
        for j = 2:length(gluedNodes{i})
            indices = (EtaMesh == gluedNodes{i}(j));
            EtaMesh(indices) = parentIdx;
        end
    end
    
    n_en = q+1;
    Fvalues = zeros(n_en,noElemsEta);
    indices = zeros(n_en,noElemsEta);
    
    [W1D,Q1D] = gaussianQuadNURBS(q+1); 
    for e = 1:noElemsEta
        idEta = indexEta(e,1);   % the index matrix is made in generateIGA3DMesh

        Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]

        J_2 = 0.5*(Eta_e(2)-Eta_e(1));

        sctrEta = xi1Nodes(EtaMesh(e,:));          %  element scatter vector
        n_en = length(sctrEta);
        pts = controlPts(sctrEta,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W1D,1)
            pt = Q1D(gp,:);
            wt = W1D(gp);

            eta  = parent2ParametricSpace(Eta_e, pt(1));
            
            [R_fun, dRdeta] = NURBS1DBasis(eta, q, Eta, weights(xi1Nodes));

            J = pts'*dRdeta';
            crossProd = J'; 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            tangent = crossProd/norm(crossProd);
            normal = [tangent(2), -tangent(1)];
            deriv = Grad(v).'*normal';

            f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrEta';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end


keyboard
if applyNeumannAtEta0
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
    Fvalues = zeros(n_en,noElemsXi);
    indices = zeros(n_en,noElemsXi);
    
    [W1D,Q1D] = gaussianQuadNURBS(p+1); 
    for e = 1:noElemsXi
        idXi = indexXi(e,1);   % the index matrix is made in generateIGA3DMesh

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]

        J_2 = 0.5*(Xi_e(2)-Xi_e(1));

        sctrXi = eta0Nodes(XiMesh(e,:));          %  element scatter vector

        n_en = length(sctrXi);
        pts = controlPts(sctrXi,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W1D,1)
            pt = Q1D(gp,:);
            wt = W1D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            
            [R_fun, dRdxi] = NURBS1DBasis(xi, p, Xi, weights(eta0Nodes));

            J = pts'*dRdxi';
            crossProd = J'; 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            tangent = -crossProd/norm(crossProd);
            normal = [tangent(2), -tangent(1)];
            deriv = Grad(v).'*normal';

            f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXi';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end


if applyNeumannAtEta1
    % Apply Neumann boundary condition at surface where eta = 1
    eta1Nodes = zeros(1,n);
    counter = 1;
    for i = 1:n
        eta1Nodes(counter) = n*(m-1) + i;
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
    Fvalues = zeros(n_en,noElemsXi);
    indices = zeros(n_en,noElemsXi);
 
    [W1D,Q1D] = gaussianQuadNURBS(p+1); 
    for e = 1:noElemsXi
        idXi = indexXi(e,1);   % the index matrix is made in generateIGA3DMesh

        Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]

        J_2 = 0.5*(Xi_e(2)-Xi_e(1));

        sctrXi = eta1Nodes(XiMesh(e,:));          %  element scatter vector

        n_en = length(sctrXi);
        pts = controlPts(sctrXi,:);
        f_e = zeros(n_en,1);
        for gp = 1:size(W1D,1)
            pt = Q1D(gp,:);
            wt = W1D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            
            [R_fun, dRdxi] = NURBS1DBasis(xi, p, Xi, weights(eta1Nodes));

            J = pts'*dRdxi';
            crossProd = J'; 


            v = R_fun*pts;
            x = v(1);
            y = v(2);
            tangent = -crossProd/norm(crossProd);
            normal = [-tangent(2), tangent(1)];
            deriv = Grad(v).'*normal';

            f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
        end    
        indices(:,e) = sctrXi';
        Fvalues(:,e) = f_e;
    end
    F = F + vectorAssembly(Fvalues,indices,noDofs);
end



