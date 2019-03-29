function p_h = calculateScatteredPressure(varCol, U, P_far, computeFarField)

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


surfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(1) == 0
        surfaceElements = [surfaceElements e];
    end
end

noGpXi = p_xi+1;
noGpEta = p_eta+1;

[W2D,Q2D] = gaussianQuadNURBS(noGpXi,noGpEta);  

if true
    p_h = zeros(size(P_far,1),1);

%     for i = 1:length(surfaceElements) %
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
            n = crossProd/J_1;

            Y = R_fun*pts;  

            p_h_gp = R_fun*U_sctr;

            if strcmp(BC, 'SHBC')
                dp_h_gp = -dp_inc(Y,n); 
            else
                dp_h_gp = dot3((J'\[dRdxi; dRdeta; dRdzeta]*U(sctr,:)).',n);
            end
            if computeFarField
                x_d_n = dot3(P_far, n)./norm2(P_far);
                x_d_y = dot3(P_far, Y.')./norm2(P_far);
                p_h = p_h - 1/(4*pi)*(1i*k* (p_h_gp.').*x_d_n + dp_h_gp).*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
            else
                xmy = -elementAddition(Y, -P_far);

                r = norm2(xmy);  
                p_h = p_h + ((p_h_gp.').*dPhi_kdny(xmy,r,n) - dp_h_gp.*Phi_k(r))* J_1 * J_2 * wt; 
            end
        end
    end
else %strcmp(formulation(3),'C') && norm(P_far(1,:)) > 1e3    
    if sum(abs(P_far(:,3))) > eps
        error('not implemented')
    end
    nurbs = varCol.nurbs;
    N = varCol.N;
    r_a = varCol.R_a;
    Upsilon = varCol.Upsilon;
	n_xi = nurbs.number(1);
	n_eta = nurbs.number(2);
    D = varCol.D;
    
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    uniqueZeta = unique(Zeta);
    noElementsXi = length(uniqueXi) - 1;
    noElementsEta = length(uniqueEta) - 1;
    noElementsZeta = length(uniqueZeta) - 1;
    R_far = norm(P_far(1,:));
    
    if false
        xi_arr = linspace(0,1,10000);
        eta = 0.5;
        p_h2 = zeros(numel(xi_arr),1);
        alpha_f_arr2 = zeros(numel(xi_arr),1);
    %     for i = 1:numel(xi_arr)
        parfor i = 1:numel(xi_arr)
            xi = xi_arr(i);
            R = NURBS3DBasis(xi, eta, 1, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);


            xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
            eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
            e = xi_idx + noElementsXi*(eta_idx-1) + noElementsXi*noElementsEta*(noElementsZeta-1);
            sctr = element(e,:);

            X = R*controlPts(sctr,:);
            [~, ~, alpha_f_arr2(i)] = evaluateProlateCoords(X(1),X(2),X(3),Upsilon);
            u = 0;
            for m = 1:N            
                phi_m = D(m,1)*exp(-1i*k*r_a)*r_a/R_far;

                u = u + phi_m*R*U(sctr + n_eta*n_xi*(m-1));    
            end
            p_h2(i) = u;
        end
    else
        x_0 = varCol.x_0;
        A_2 = varCol.A_2;
        xi_arr = [0, 0.5];
        eta_arr = fliplr(linspace(0,1,10000));
        p_h2 = zeros(numel(eta_arr),2);
        alpha_f_arr2 = zeros(numel(eta_arr),2);
        for j = 1:2
            xi = xi_arr(j);
        %     for i = 1:numel(xi_arr)
            parfor i = 1:numel(eta_arr)
                eta = eta_arr(i);
                R = NURBS3DBasis(xi, eta, 1, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);

                xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
                eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
                e = xi_idx + noElementsXi*(eta_idx-1) + noElementsXi*noElementsEta*(noElementsZeta-1);
                sctr = element(e,:);

                X = R*controlPts(sctr,:);
                Xt = A_2*(X-x_0)';
                [~, alpha_f_arr2(i,j), ~] = evaluateProlateCoords(Xt(1),Xt(2),Xt(3),Upsilon);
                u = 0;
                for m = 1:N            
                    phi_m = D(m,1)*exp(-1i*k*r_a)*r_a/R_far;

                    u = u + phi_m*R*U(sctr + n_eta*n_xi*(m-1));    
                end
                p_h2(i,j) = u;
            end
        end
        p_h2 = p_h2(:);
        
        alpha_f_arr2(:,2) = alpha_f_arr2(:,2) + pi;
        alpha_f_arr2 = alpha_f_arr2(:);
        [alpha_f_arr2, I] = unique(alpha_f_arr2);
        p_h2 = p_h2(I);
    end
    
    alpha_f_arr2(1) = 0;
    alpha_f_arr2(end) = 2*pi;
    
    alpha_f_arr = atan2(P_far(:,2),P_far(:,1));
    alpha_f_arr(alpha_f_arr < 0) = alpha_f_arr(alpha_f_arr < 0) + 2*pi;
    
    p_h = interp1(alpha_f_arr2,p_h2,alpha_f_arr,'spline');
end



