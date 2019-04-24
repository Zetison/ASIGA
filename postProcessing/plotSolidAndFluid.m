
Xi = solid.knots{1};
Eta = solid.knots{2};
Zeta = solid.knots{3};
noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));

extraXiPts = floor(40/(noUniqueXiKnots-1)); % .. per element
extraEtaPts = floor(100/(noUniqueEtaKnots-1)); % .. per element
extraZetaPts = floor(0/(noUniqueZetaKnots-1)); % .. per element

vtfFileName = ['../graphics/GLview/scatteringOnElasticSphere/' num2str(sample) '_FEM' num2str(analyzeFEMcase)];
options = {'name',vtfFileName, ...
             'plotDisplacementVectors',1,...
             'plotXdisplacement',0,...
             'plotYdisplacement',0,...
             'plotZdisplacement',0,...
             'plotSphericalRadialDisplacement',0,...
             'plotVonMisesStress',1,...
             'plotTimeOscillation', 1,...
             'plotMagnitudeDisplacement',0,...
             'addDynamicScalars', 1};


maxJacobian = 1700;

weights = weights_solid;

[nodes_solid, noNodes_solid, visElements_solid, cornerNode_solid, ...
 noXiKnots_solid, noEtaKnots_solid, noZetaKnots_solid] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems_solid, solid);

stress = zeros(noNodes_solid,6);
displacement   = zeros(noNodes_solid,3);
sigma_v = zeros(noNodes_solid,1);

for e = 1:noElems_solid
    idXi = index_solid(e,1);
    idEta = index_solid(e,2);
    idZeta = index_solid(e,3);
    
    % Modify element range such that the knotspan index do not change over
    % an element (with value 1e-12)
    
    Xi_e = elRangeXi_solid(idXi,:); % [xi_i,xi_i+1]
    Xi_eV = linspace(Xi_e(1),Xi_e(end)-1e-12,2+extraXiPts);
    
    Eta_e = elRangeEta_solid(idEta,:); % [eta_j,eta_j+1]
    Eta_eV = linspace(Eta_e(1),Eta_e(end)-1e-12,2+extraEtaPts);
    
    Zeta_e = elRangeZeta_solid(idZeta,:); % [zeta_k,zeta_k+1]
    Zeta_eV = linspace(Zeta_e(1),Zeta_e(end)-1e-12,2+extraZetaPts);
    
    sctr = element_solid(e,:);          %  element scatter vector
    sctrB = [sctr sctr+noCtrlPts sctr+2*noCtrlPts]; % scatters a B matrix
    n_en = length(sctr);
    
    B = zeros(6,3*n_en);
    pts = controlPts_solid(sctr,:);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
        
    for idxZeta=1:2+extraZetaPts  
        for idxEta=1:2+extraEtaPts
            for idxXi = 1:2+extraXiPts       
                xi = Xi_eV(idxXi);       
                eta = Eta_eV(idxEta);    
                zeta = Zeta_eV(idxZeta);    
                nodeID = cornerNode_solid(e) + idxXi-1   +   (idxEta-1)*noXiKnots_solid + (idxZeta-1)*noEtaKnots_solid*noXiKnots_solid;                                      
                
                [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_solid, q_solid, r_solid, Xi, Eta, Zeta, weights_solid);
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
                
                % B matrix
                dRdx = dRdX(1,:)*elemDisp;
                dRdy = dRdX(2,:)*elemDisp;
                dRdz = dRdX(3,:)*elemDisp;
                
                strain = calculateStrainVector(dRdx, dRdy, dRdz);
                
                
                
                stress(nodeID,:) = C*strain;
                
            end
        end
    end
end

data.nodes = nodes_solid;
data.visElements = visElements_solid;
data.displacement = displacement;
data.stress = stress;
if exist('omega','var')
    data.omega = omega;
else
    data.omega = 1;
end


makeVTFfile_new(data, options);




vtfFileName = ['../graphics/GLview/scatteringOnElasticSphere/' num2str(sample) 'fluid'];

options = {'name',vtfFileName, ...
             'plotScalarField',1};

activateFluidParameters

noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));

extraXiPts = floor(40/(noUniqueXiKnots-1)); % .. per element
extraEtaPts = floor(40/(noUniqueEtaKnots-1)); % .. per element
extraZetaPts = floor(0/(noUniqueZetaKnots-1)); % .. per element
plotScalars_new