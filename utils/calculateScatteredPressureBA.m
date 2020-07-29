function p_h = calculateScatteredPressureBA(varCol, U, P_far, useExtraQuadPts, computeFarField)
error('Depricated, use calculateScattere instead')
formulation = varCol.formulation;
solveForPtot = varCol.solveForPtot;
SHBC = strcmp(varCol.BC, 'SHBC');
if SHBC
    if solveForPtot
        homNeumanCond = true;
        dpdn = @(x,n) 0;
    else
        homNeumanCond = false;
        dp_inc = varCol.dp_inc;
        dpdn = @(x,n) -dp_inc(x,n);
    end
else
    homNeumanCond = false;
    dpdn = varCol.dpdn;
end
degree = varCol.degree; % assume degree is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
d_p = varCol.patches{1}.nurbs.d_p;
extraGP = varCol{1}.extraGP;
k = varCol.k;
if useExtraQuadPts
    noGpXi = degree+1+5;
    noGpEta = degree+1+5;
else
    noGpXi = degree+1+extraGP;
    noGpEta = degree+1+extraGP;
end
[Q, W] = gaussTensorQuad([noGpXi,noGpEta]); 

switch formulation
    case 'SL2E'
                if strcmp(varCol.coreMethod, 'XI')
                    useEnrichedBfuns = true;
                    d_vec = varCol.d_vec;
                else
                    useEnrichedBfuns = false;
                    d_vec = NaN;
                end

                % [W2D,Q2D] = gaussianQuadNURBS(20,20);

                p_h = zeros(size(P_far,1),1);
        %         for e = 1:noElems %
        parfor e = 1:noElems
            patch = pIndex(e);
            knots = knotVecs{patch};
            Xi_e = zeros(d_p,2);
            for ii = 1:d_p
                Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
            end

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:);
            J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;


            xi = parent2ParametricSpace(Xi_e, Q);
            I = findKnotSpans(degree, xi(1,:), knots);
            R = NURBSbasis(I, xi, degree, knots, wgts);
            [J_1,crossProd] = getJacobian(R,pts,d_p);

            for gp = 1:size(W,1)
                wt = W(gp);
                n = crossProd.'/J_1;

                Y = R{1}(qp,:)*pts;

                xmy = -elementAddition(Y, -P_far);

                r = norm2(xmy);    
                if useEnrichedBfuns
                    R{1} = R{1}*exp(1i*k*dot(d_vec, Y));
                end

                p_h_gp = R{1}(qp,:)*U(sctr,:);
                if computeFarField
                    x_d_n = dot3(P_far, n.')./norm2(P_far);
                    x_d_y = dot3(P_far, Y.')./norm2(P_far);
                    p_h = p_h + -1/(4*pi)*1i*k* (p_h_gp.').*x_d_n.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                else
                    p_h = p_h + (p_h_gp.').*dPhi_kdny(xmy,r,n,k)* J_1 * J_2 * wt;  
                end
                if ~homNeumanCond
                    dp_h_gp = dpdn(Y, n).';
                    if computeFarField
                        p_h = p_h + -1/(4*pi)*dp_h_gp.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                    else
                        p_h = p_h + -dp_h_gp.*Phi_k(r,k)* J_1 * J_2 * wt;  
                    end
                end
            end
        end
    case 'VL2E'

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
            for gp = 1:size(W,1)
                pt = Q(gp,:);
                wt = W(gp);

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

                if computeFarField
                    x_d_n = dot3(P_far, n.')./norm2(P_far);
                    x_d_y = dot3(P_far, Y.')./norm2(P_far);
                    p_h = p_h + -1/(4*pi)*1i*k* (p_h_gp.').*x_d_n.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                else
                    xmy = -elementAddition(Y, -P_far);
                    r = norm2(xmy);    
                    p_h = p_h + (p_h_gp.').*dPhi_kdny(xmy,r,n,k)* J_1 * J_2 * wt;  
                end
                if ~homNeumanCond
                    if SHBC
                        dp_h_gp = dpdn(Y, n);
                    else
                        dp_h_gp = dot3((J'\[dRdxi; dRdeta; dRdzeta]*U(sctr,:)).',n);
                    end
                    if computeFarField
                        p_h = p_h + -1/(4*pi)*dp_h_gp.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                    else
                        p_h = p_h + -dp_h_gp.*Phi_k(r,k)* J_1 * J_2 * wt;  
                    end
                end
            end
        end
end





