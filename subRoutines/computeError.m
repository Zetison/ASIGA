energyError_e = 0;
normOfAnalytic = 0;

noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));

% vtfFileName = '../graphics/GLview/pinchedHemisphere';
extraXiPts = 2^(noAnalysisArr(degreeElev+1) - M) - 1; % .. per element
extraEtaPts = 2^(noAnalysisArr(degreeElev+1) - M) - 1; % .. per element
extraZetaPts = 2^(noAnalysisArr(degreeElev+1) - M) - 1; % .. per element

[nodes, noNodes, visElements, cornerNode, ...
          noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, solid);

[W3D,Q3D] = gaussianQuadNURBS(p+3,q+3,r+3); 

for e = 1:noElems
    idXi    = index(e,1);
    idEta    = index(e,2);
    idZeta    = index(e,3);
    
    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Xi_eV = linspace(Xi_e(1),Xi_e(end),2+extraXiPts);
    
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Eta_eV = linspace(Eta_e(1),Eta_e(end),2+extraEtaPts);
    
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    Zeta_eV = linspace(Zeta_e(1),Zeta_e(end),2+extraZetaPts);
    
    sctr   = element(e,:);          %  element scatter vector
        
    J_2 = 0.125*(Xi_eV(2)-Xi_eV(1))*(Eta_eV(2)-Eta_eV(1))*(Zeta_eV(2)-Zeta_eV(1));
    
    pts    = controlPts(sctr,:);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
    
    % loop over Gauss points
    for idxZeta=1:2+extraZetaPts-1
        Zeta_eVe = [Zeta_eV(idxZeta) Zeta_eV(idxZeta+1)];
        for idxEta=1:2+extraEtaPts-1
            Eta_eVe = [Eta_eV(idxEta) Eta_eV(idxEta+1)];
            for idxXi = 1:2+extraXiPts-1   
                Xi_eVe = [Xi_eV(idxXi) Xi_eV(idxXi+1)];     
                for gp = 1:size(W3D,1)
                    pt = Q3D(gp,:);
                    wt = W3D(gp);

                    % compute coords in parameter space
                    xi = parent2ParametricSpace(Xi_eVe,  pt(1));
                    eta = parent2ParametricSpace(Eta_eVe, pt(2));
                    zeta = parent2ParametricSpace(Zeta_eVe,pt(3));

                    % compute derivative of basis functions w.r.t parameter coord
                    [~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);

                    J = pts'*[dRdxi' dRdeta' dRdzeta'];
                    J_1 = det(J);        

                    dRdX = J'\[dRdxi; dRdeta; dRdzeta];

                    % Compute error contribution
                    dRdx = dRdX(1,:)*elemDisp;
                    dRdy = dRdX(2,:)*elemDisp;
                    dRdz = dRdX(3,:)*elemDisp;

                    strain_h = calculateStrainVector(dRdx, dRdy, dRdz);

                    v = evaluateNURBS(solid, [xi, eta, zeta]);
                    x = v(1);
                    y = v(2);
                    z = v(3);
                    theta = atan2(y,x);
                    radius = sqrt(x^2+y^2);

                    strain_a = strain_u(radius,theta,z);
                    strain_e = strain_a - strain_h;


                    energyError_e  = energyError_e  + strain_e' * C * strain_e * abs(J_1) * J_2 * wt;
                    normOfAnalytic = normOfAnalytic + strain_a' * C * strain_a * abs(J_1) * J_2 * wt;
                end
            end
        end
    end
end
normOfAnalytic
% with c_i = 1 i=1,2,3 normOfAnalytic = 3.593981996e-11

relativeEnergyNorm = sqrt(energyError_e)/sqrt(normOfAnalytic);



%% The following will compute the max norm
% noPtsXi = 5;
% noPtsEta = noPtsXi;
% notPtsZeta = noPtsXi;
% 
% 
% xiArray   = linspace(0,	1, noPtsXi);
% etaArray  = linspace(0,	1, noPtsEta);
% zetaArray = linspace(0, 1, notPtsZeta);
% 
% maxNorm = zeros(noPtsXi*noPtsEta*notPtsZeta,1);
% 
% counter = 1;
% for k = 1:noPtsXi
%     for j = 1:noPtsEta
%         for i = 1:notPtsZeta
%             v = evaluateNURB(xiArray(i), etaArray(j), zetaArray(k), solid);
% 
%             uAtPoint = u_analytic(v(1),v(2),v(3));
%             maxNorm(counter) = norm(uAtPoint - numericalSolEval(xiArray(i), etaArray(j), zetaArray(k), p, q, r, Xi, Eta, Zeta, weights, Ux, Uy, Uz),inf);
%             counter = counter + 1;
%         end
%     end
% end
% maxNormArray(analysisCounter) = max(maxNorm);

