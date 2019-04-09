function p_h = calculateScatteredPressureBA(varCol, U, P_far, useExtraQuadPts, computeFarField)

switch varCol.formulation
    case 'SL2E'
        p_xi = varCol.degree(1); % assume p_xi is equal in all patches
        p_eta = varCol.degree(2); % assume p_eta is equal in all patches

        index = varCol.index;
        noElems = varCol.noElems;
        elRangeXi = varCol.elRange{1};
        elRangeEta = varCol.elRange{2};
        element = varCol.element;
        element2 = varCol.element2;
        weights = varCol.weights;
        controlPts = varCol.controlPts;
        knotVecs = varCol.knotVecs;
        pIndex = varCol.pIndex;
        
        k = varCol.k;
        BC = varCol.BC;
        SHBC = strcmp(BC, 'SHBC');
        if SHBC
            dpdn = 0;
            dp_inc = varCol.dp_inc;
        else
            dpdn = varCol.dpdn;
            dp_inc = NaN;
        end

        Phi_k = varCol.Phi_k;
        dPhi_kdny = varCol.dPhi_kdny;


        if strcmp(varCol.coreMethod, 'XI')
            useEnrichedBfuns = true;
            d_vec = varCol.d_vec;
        else
            useEnrichedBfuns = false;
            d_vec = NaN;
        end
        if useExtraQuadPts
            noGpXi = p_xi+5;
            noGpEta = p_eta+5;
        else
            noGpXi = p_xi+1;
            noGpEta = p_eta+1;
        end

        [W2D,Q2D] = gaussianQuadNURBS(noGpXi,noGpEta);  
        % [W2D,Q2D] = gaussianQuadNURBS(20,20);

        p_h = zeros(size(P_far,1),1);
%         for e = 1:noElems %
        parfor e = 1:noElems
            patch = pIndex(e); % New
            Xi = knotVecs{patch}{1}; % New
            Eta = knotVecs{patch}{2}; % New

            idXi = index(e,1);
            idEta = index(e,2);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));


            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);
                wt = W2D(gp);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));

                [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd.'/J_1;

                Y = R*pts;

                xmy = -elementAddition(Y, -P_far);

                r = norm2(xmy);    
                if useEnrichedBfuns
                    R = R*exp(1i*k*dot(d_vec, Y));
                end

                p_h_gp = R*U(sctr,:);
                if SHBC
                    dp_h_gp = -dp_inc(Y,n);
                else
                    dp_h_gp = dpdn(Y,n);   
                end
                if computeFarField
                    x_d_n = dot3(P_far, n.')./norm2(P_far);
                    x_d_y = dot3(P_far, Y.')./norm2(P_far);
                    p_h = p_h - 1/(4*pi)*(1i*k* (p_h_gp.').*x_d_n + dp_h_gp).*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                else
                    p_h = p_h + ((p_h_gp.').*dPhi_kdny(xmy,r,n) - dp_h_gp.*Phi_k(r))* J_1 * J_2 * wt; 
                end
            end
        end
    case 'VL2E'
        p_xi = varCol.degree(1); % assume p_xi is equal in all patches
        p_eta = varCol.degree(2); % assume p_eta is equal in all patches
        p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches

        index = varCol.index;
        noElems = varCol.noElems;
        elRangeXi = varCol.elRange{1};
        elRangeEta = varCol.elRange{2};
        elRangeZeta = varCol.elRange{3};
        element = varCol.element;
        element2 = varCol.element2;
        weights = varCol.weights;
        controlPts = varCol.controlPts;
        knotVecs = varCol.knotVecs;
        pIndex = varCol.pIndex;
        
        BC = varCol.BC;
        k = varCol.k;

        dp_inc = varCol.dp_inc;

        Phi_k = varCol.Phi_k;
        dPhi_kdny = varCol.dPhi_kdny;


        % Find elements on the inner surface for evaluation of
        % backscattered pressure in far field
        surfaceElements = [];
        for e = 1:noElems
            idZeta = index(e,3);
            Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]                    
            if Zeta_e(1) == 0
                surfaceElements = [surfaceElements e];
            end
        end
        if useExtraQuadPts
            noGpXi = p_xi+5;
            noGpEta = p_eta+5;
        else
            noGpXi = p_xi+1;
            noGpEta = p_eta+1;
        end

        [W2D,Q2D] = gaussianQuadNURBS(noGpXi,noGpEta);  


        p_h = zeros(size(P_far,1),1);

%         for i = 1:length(surfaceElements) %
        parfor i = 1:length(surfaceElements) %[8 7 3 4 5 6 2 1]% 

            e = surfaceElements(i);
            patch = pIndex(e); % New
            Xi = knotVecs{patch}{1}; % New
            Eta = knotVecs{patch}{2}; % New
            Zeta = knotVecs{patch}{3}; % New
            idXi = index(e,1);
            idEta = index(e,2);

            Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
            Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            sctr = element(e,:);
            U_sctr = U(sctr,:);

            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New
            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);
                wt = W2D(gp);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));
                zeta = 0;

                [R_fun, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

                J = pts'*[dRdxi' dRdeta' dRdzeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd.'/J_1;

                Y = R_fun*pts;  

                p_h_gp = R_fun*U_sctr;

                if strcmp(BC, 'SHBC')
                    dp_h_gp = -dp_inc(Y,n); 
                else
                    dp_h_gp = dot3((J'\[dRdxi; dRdeta; dRdzeta]*U(sctr,:)).',n);
                end
                if computeFarField
                    x_d_n = dot3(P_far, n.')./norm2(P_far);
                    x_d_y = dot3(P_far, Y.')./norm2(P_far);
                    p_h = p_h - 1/(4*pi)*(1i*k* (p_h_gp.').*x_d_n + dp_h_gp).*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                else
                    xmy = -elementAddition(Y, -P_far);

                    r = norm2(xmy);  
                    p_h = p_h + ((p_h_gp.').*dPhi_kdny(xmy,r,n) - dp_h_gp.*Phi_k(r))* J_1 * J_2 * wt; 
                end
            end
        end
end





