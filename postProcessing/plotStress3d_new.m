
% build visualization B8 mesh
Xi = nurbs.knots{1};
Eta = nurbs.knots{2};
Zeta = nurbs.knots{3};

[nodes, noNodes, visElements, cornerNode, ...
 noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, nurbs);

stress = zeros(noNodes,6);
displacement   = zeros(noNodes,3);
orgJacobi = zeros(noNodes,1);

noElems = (noUniqueXiKnots-1)*(noUniqueEtaKnots-1)*(noUniqueZetaKnots-1);
for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    % Modify element range such that the knotspan index do not change over
    % an element (with value 1e-12)
    
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Xi_eV = linspace(Xi_e(1),Xi_e(end)-1e-12,2+extraXiPts);
    
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Eta_eV = linspace(Eta_e(1),Eta_e(end)-1e-12,2+extraEtaPts);
    
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    Zeta_eV = linspace(Zeta_e(1),Zeta_e(end)-1e-12,2+extraZetaPts);
    
    sctr = element(e,:);          %  element scatter vector
    sctrB = [sctr sctr+noCtrlPts sctr+2*noCtrlPts]; % scatters a B matrix
    n_en = length(sctr);
    
    B = zeros(6,3*n_en);
    pts = controlPts(sctr,:);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
        
    for idxZeta=1:2+extraZetaPts  
        for idxEta=1:2+extraEtaPts
            for idxXi = 1:2+extraXiPts           
                xi = Xi_eV(idxXi);   
                eta = Eta_eV(idxEta); 
                zeta = Zeta_eV(idxZeta);                                             
                
                nodeID = cornerNode(e) + idxXi-1   +   (idxEta-1)*noXiKnots + (idxZeta-1)*noEtaKnots*noXiKnots;
                
                [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                displacement(nodeID,:) = R_fun*elemDisp;   
                % compute the jacobian of physical and parameter domain mapping
                % then the derivative w.r.t spatial physical coordinates
                
                J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                if cond(J) == inf
                    orgJacobi(nodeID) = 6e11; % -(cond(J) - maxJacobian);
                else
                    orgJacobi(nodeID) = cond(J);  % -(cond(J) - maxJacobian);
                end
                    %                 
%                 if idxXi == ceil((2+extraXiPts)/2) && idxEta == ceil((2+extraEtaPts)/2) && idxZeta == ceil((2+extraZetaPts)/2)
%                     keyboard
%                 end
                if cond(J) > maxJacobian
                    if idxXi == 1
                        xi = Xi_eV(2);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    elseif idxXi == 2+extraXiPts
                        xi = Xi_eV(idxXi-1);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    end
                end
                if cond(J) > maxJacobian
                    if idxEta == 1
                        eta = Eta_eV(2);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    elseif idxEta == 2+extraEtaPts
                        eta = Eta_eV(idxEta-1);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    end
                end
                if cond(J) > maxJacobian
                    if idxZeta == 1
                        zeta = Zeta_eV(2);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    elseif idxZeta == 2+extraZetaPts
                        zeta = Zeta_eV(idxZeta-1);
                        [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    end
                end
                
                % Jacobian inverse and spatial derivatives                
                dRdX = J'\[dRdxi; dRdeta; dRdzeta];
                
                
                dRdx = dRdX(1,:)*elemDisp;
                dRdy = dRdX(2,:)*elemDisp;
                dRdz = dRdX(3,:)*elemDisp;
                
                strain = calculateStrainVector(dRdx, dRdy, dRdz);
                
                
                
                stress(nodeID,:) = C*strain;
            end
        end
    end
end

data.nodes = nodes;
data.visElements = visElements;
data.displacement = displacement;
data.stress = stress;
data.orgJacobi = orgJacobi;
if exist('omega','var')
    data.omega = omega;
else
    data.omega = 1;
end

makeVTKfile_new(data, options);



