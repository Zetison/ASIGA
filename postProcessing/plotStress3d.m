
% build visualization B8 mesh

buildVisualization3dMesh

stress = zeros(noElems,size(elementV,2),6);
disp   = zeros(noElems,size(elementV,2),3);

sigma_v = zeros(noElems,8);

for e=1:noElems
    idXi    = index(e,1);
    idEta    = index(e,2);
    idZeta    = index(e,3);
    
    Xi_e    = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e   = elRangeEta(idEta,:); % [eta_j,eta_j+1]
    Zeta_e  = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
    
    % Modify element range such that the knotspan index do not change over
    % an element
    
    sctr   = element(e,:);          %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts sctr+2*noCtrPts]; % scatters a B matrix
    n_en     = length(sctr);
    
    B      = zeros(6,3*n_en);
    pts    = controlPts(sctr,:);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
    
    % loop over Gauss points
    
    gp = 1;
    for idxXi=1:2
        zeta = Zeta_e(idxXi);
        for idxEta=1:2
            eta  = Eta_e(idxEta);
            for idxZeta=1:2                
                xi   = Xi_e(idxZeta);                                                
                
                [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
                % compute the jacobian of physical and parameter domain mapping
                % then the derivative w.r.t spatial physical coordinates
                
                J  = pts' * [dRdxi' dRdeta' dRdzeta'];
                
                % Jacobian inverse and spatial derivatives
                
                dRdX = J'\[dRdxi; dRdeta; dRdzeta];
                
                % B matrix
                dRdx = dRdX(1,:)*elemDisp;
                dRdy = dRdX(2,:)*elemDisp;
                dRdz = dRdX(3,:)*elemDisp;
                
                strain = calculateStrainVector(dRdx, dRdy, dRdz);
                
                stress(e,gp,:)  = C*strain;
                disp  (e,gp,:)  = R_fun*elemDisp;   
                
                
                gp = gp + 1;
            end
        end
        sigma_v(e,:) = sqrt((stress(e,:,1)-stress(e,:,2)).^2 + (stress(e,:,2)-stress(e,:,3)).^2 + (stress(e,:,1)-stress(e,:,3)).^2 ...
            + 6*(stress(e,:,4).^2 + stress(e,:,5).^2 + stress(e,:,6).^2)/2);
    end
    
    % disp stored in IGA element connectivity
    % change positions according to standard FE connectivity
    
    col3 = disp(e,3,:);
    col4 = disp(e,4,:);
    col7 = disp(e,7,:);
    col8 = disp(e,8,:);
    
    disp(e,3,:) = col4;
    disp(e,4,:) = col3;
    disp(e,7,:) = col8;
    disp(e,8,:) = col7;
    
    col3 = stress(e,3,:);
    col4 = stress(e,4,:);
    col7 = stress(e,7,:);
    col8 = stress(e,8,:);
    
    stress(e,3,:) = col4;
    stress(e,4,:) = col3;
    stress(e,7,:) = col8;
    stress(e,8,:) = col7;
    
    col3 = sigma_v(e,3);
    col4 = sigma_v(e,4);
    col7 = sigma_v(e,7);
    col8 = sigma_v(e,8);
    
    sigma_v(e,3,:) = col4;
    sigma_v(e,4,:) = col3;
    sigma_v(e,7,:) = col8;
    sigma_v(e,8,:) = col7;
end
% displacements
noNodes = size(nodes,1);

dispX   = zeros(noNodes,1);
dispY   = zeros(noNodes,1);
dispZ   = zeros(noNodes,1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:8
        nid = connect(in);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
        dispZ(nid) = disp(e,in,3);
    end
end

makeVTFfile

