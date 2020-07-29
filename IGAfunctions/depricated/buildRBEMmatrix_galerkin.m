function [A, FF] = buildRBEMmatrix_galerkin(varCol)
error('Depricated. Use buildGBEMmatrix instead')
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
noElemsPatch = varCol.noElemsPatch;
noPatches = varCol.noPatches;
patches = varCol.patches;
extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;

dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;



model = varCol.model;
k = varCol.k;
alpha = 1i/k;

switch varCol.formulation(3:end)
    case 'BM'
        useCBIE = true;
        useHBIE = true;
    case 'CBIE'
        useCBIE = true;
        useHBIE = false;
    case 'HBIE'
        useCBIE = false;
        useHBIE = true;
    otherwise
        error('Formulation not implemented')
end


Phi_0 = @(r)           1/(4*pi*r);
Phi_k = @(r) exp(1i*k*r)/(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)/r^2*             (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)/r^2*(1 - 1i*k*r)*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

no_angles = length(varCol.alpha_s);
p_inc = varCol.p_inc;
dp_inc = varCol.dp_inc;

exteriorProblem = true;
if exteriorProblem
    sgn = -1;
else
    sgn = 1;
end


n_en = (p_xi+1)*(p_eta+1);
% 
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);
p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
% [W2D,Q2D] = gaussianQuadNURBS(2*p_xi+1,2*p_eta+1);
% p_max = max(p_xi,p_eta);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(6*p_max+1,6*p_max+1);

idxRow = zeros(n_en, noElems);
Avalues = zeros(n_en, noDofs, noElems); 
Fvalues = zeros(n_en, noElems, no_angles); 
parfor e_x = 1:noElems
% for e_x = 1:noElems
    patch_x = pIndex(e_x); % New
    Xi_x = knotVecs{patch_x}{1}; % New
    Eta_x = knotVecs{patch_x}{2}; % New

    idXi_x = index(e_x,1);
    idEta_x = index(e_x,2);

    Xi_e_x = elRangeXi(idXi_x,:);
    Eta_e_x = elRangeEta(idEta_x,:);

    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts_x = weights(element2(e_x,:)); % New 

    F_e = zeros(n_en, no_angles);
    A_e_temp = zeros(n_en, n_en, noElems);
    
    J_2 = 0.25*(Xi_e_x(2)-Xi_e_x(1))*(Eta_e_x(2)-Eta_e_x(1));
    
    for gp_x = 1:size(W2D,1)
        pt_x = Q2D(gp_x,:);
        wt_x = W2D(gp_x);

        xi_x  = parent2ParametricSpace(Xi_e_x, pt_x(1));
        eta_x = parent2ParametricSpace(Eta_e_x,pt_x(2));
        [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi_x, Eta_x, wgts_x);

        J = [dR_xdxi; dR_xdeta]*pts_x;
        m_1 = J(1,:);
        m_2 = J(2,:);
        crossProd = cross(m_1,m_2);
        J_1 = norm(crossProd);

        x = R_x*pts_x;
        fact_x = J_1*J_2*wt_x;

        h_xi = norm(m_1);
        h_eta = norm(m_2);
        e_xi = m_1/h_xi;
        e_eta = m_2/h_eta;

        v_1 = e_xi;
        nx = crossProd'/J_1;
        v_3 = nx.';
        v_2 = cross(v_3,v_1);
        cosT = dot(e_xi,e_eta);
        sinT = dot(v_2,e_eta);
        dXIdv = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
        
        if abs(nx(1)) < 1/sqrt(2)
            d1 = sqrt(3)/2*cross([1;0;0],nx)/sqrt(1-nx(1)^2) - nx/2;
        else
            d1 = sqrt(3)/2*cross([0;1;0],nx)/sqrt(1-nx(2)^2) - nx/2;
        end
        d2 = d1+nx;
        Psi3 = 0;
        Psi4 = 0;
        dPsi3dny = 0;
        dPsi4dny = 0;
        
        Psi1_integral = 0;
        Psi2_integral = 0;
        Psi3_integral = 0;
        Psi4_integral = 0;
        dPsi1dny_integral = 0;
        dPsi2dny_integral = 0;
        
        dPhi_0dny_integral = 0;
        d2Phi_0dnxdny_integral = 0;
        ugly_integral = zeros(3,1);
        idxCol = zeros(n_en, noElems);
    
        for e_y = 1:noElems   
            patch_y = pIndex(e_y); % New
            Xi_y = knotVecs{patch_y}{1}; % New
            Eta_y = knotVecs{patch_y}{2}; % New

            idXi_y = index(e_y,1);
            idEta_y = index(e_y,2);

            Xi_e_y = elRangeXi(idXi_y,:);
            Eta_e_y = elRangeEta(idEta_y,:);

            sctr_y = element(e_y,:);
            pts_y = controlPts(sctr_y,:);
            wgts_y = weights(element2(e_y,:)); % New 

            if useCBIE
                CBIE = zeros(1, n_en);
            end
            if useHBIE
                HBIE = zeros(1, n_en);
            end
            if e_x == e_y
                xi_x_t = parametric2parentSpace(Xi_e_y, xi_x);
                eta_x_t = parametric2parentSpace(Eta_e_y, eta_x);
                theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
                theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
                theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
                theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

                J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));

                for area = {'South', 'East', 'North', 'West'}
                    switch area{1}
                        case 'South'
                            thetaRange = [theta_x3 theta_x4];
                        case 'East'
                            thetaRange = [theta_x4 theta_x1];
                        case 'North'
                            thetaRange = [theta_x1 theta_x2];
                        case 'West'
                            if theta_x3 < 0
                                thetaRange = [theta_x2 theta_x3+2*pi];
                            else
                                thetaRange = [theta_x2 theta_x3];
                            end
                    end
                    for gp_y = 1:size(W2D_2,1)
                        pt_y = Q2D_2(gp_y,:);
                        wt_y = W2D_2(gp_y);

                        rho_t = parent2ParametricSpace([0, 1],   pt_y(1));
                        theta = parent2ParametricSpace(thetaRange,pt_y(2));
                        switch area{1}
                            case 'South'
                                rho_hat = (-1 - eta_x_t)/sin(theta);
                            case 'East'
                                rho_hat = ( 1 - xi_x_t)/cos(theta);
                            case 'North'
                                rho_hat = ( 1 - eta_x_t)/sin(theta);
                            case 'West'
                                rho_hat = (-1 - xi_x_t)/cos(theta);
                        end
                        rho = rho_hat*rho_t;

                        xi_y_t  = xi_x_t + rho*cos(theta);
                        eta_y_t = eta_x_t + rho*sin(theta);

                        J_3_y = rho;
                        J_4_y = rho_hat;
                        J_5_y = 0.25*(thetaRange(2)-thetaRange(1));

                        xi_y = parent2ParametricSpace(Xi_e_y, xi_y_t);
                        eta_y = parent2ParametricSpace(Eta_e_y, eta_y_t);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                        J_y = [dR_ydxi; dR_ydeta]*pts_y;
                        crossProd_y = cross(J_y(1,:),J_y(2,:));
                        J_1_y = norm(crossProd_y);
                        ny = crossProd_y'/J_1_y;

                        y = R_y*pts_y;
                        
                        xmy = x-y;
                        r = norm(xmy);
                        fact_y = J_1_y*J_2_y*J_3_y*J_4_y*J_5_y*wt_y;
                        if useCBIE                            
                            exp1 = exp(1i*k*dot(d1,y-x));
                            exp2 = exp(1i*k*dot(d2,y-x));
                            Psi2 = 1i*(exp1-exp2)/k;
                            Psi1 = (exp1+exp2)/2;
                            dPsi2dny = dot(d2,ny)*exp2 - dot(d1,ny)*exp1;
                            dPsi1dny = 1i*k*(dot(d1,ny)*exp1+dot(d2,ny)*exp2)/2;
                            
                            dPhi_kTemp = dPhi_kdny(xmy,r,ny);
                            Phi_kTemp = Phi_k(r);

                            Psi1_integral     = Psi1_integral    + Psi1*dPhi_kTemp*fact_y; 
                            Psi2_integral     = Psi2_integral    + Psi2*dPhi_kTemp*fact_y; 
                            Psi3_integral     = Psi3_integral    + (dPsi3dny*Phi_kTemp-Psi3*dPhi_kTemp)*fact_y; 
                            Psi4_integral     = Psi4_integral    + (dPsi4dny*Phi_kTemp-Psi4*dPhi_kTemp)*fact_y; 
                            dPsi1dny_integral = dPsi1dny_integral + dPsi1dny*Phi_kTemp*fact_y;
                            dPsi2dny_integral = dPsi2dny_integral + dPsi2dny*Phi_kTemp*fact_y;

                            CBIE = CBIE + dPhi_kTemp*R_y*fact_y;
                        end
                    end  
                end
            else
                nurbs = patches{patch_y}.nurbs;
                x_1 = evaluateNURBS(nurbs, [Xi_e_y(1)+eps, Eta_e_y(1)+eps]).';
                x_2 = evaluateNURBS(nurbs, [Xi_e_y(2)-eps, Eta_e_y(1)+eps]).';
                x_3 = evaluateNURBS(nurbs, [Xi_e_y(2)-eps, Eta_e_y(2)-eps]).';
                x_4 = evaluateNURBS(nurbs, [Xi_e_y(1)+eps, Eta_e_y(2)-eps]).';
                x_5 = evaluateNURBS(nurbs, [mean(Xi_e_y),  mean(Eta_e_y)]).';

                l = norm(x-x_5);
                h_1 = norm(x_1-x_3);
                h_2 = norm(x_2-x_4);
                h = max(h_1,h_2);
                n_div = round(3*h/l + 1);
                
                Xi_e_arr  = linspace(Xi_e_y(1),Xi_e_y(2),n_div+1);
                Eta_e_arr = linspace(Eta_e_y(1),Eta_e_y(2),n_div+1);
                for i_eta = 1:n_div
                    Eta_e_sub = Eta_e_arr(i_eta:i_eta+1);
                    for i_xi = 1:n_div
                        Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);
                        J_2_y = 0.25*(Xi_e_sub(2)-Xi_e_sub(1))*(Eta_e_sub(2)-Eta_e_sub(1));
                        for gp_y = 1:size(W2D,1)
                            pt_y = Q2D(gp_y,:);
                            wt_y = W2D(gp_y);

                            xi_y  = parent2ParametricSpace(Xi_e_sub, pt_y(1));
                            eta_y = parent2ParametricSpace(Eta_e_sub,pt_y(2));
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                            J_y = [dR_ydxi; dR_ydeta]*pts_y;
                            crossProd_y = cross(J_y(1,:),J_y(2,:));
                            J_1_y = norm(crossProd_y);
                            ny = crossProd_y'/J_1_y;

                            y = R_y*pts_y;

                            xmy = x-y;
                            r = norm(xmy);
                            fact_y = J_1_y*J_2_y*wt_y;
                            if useCBIE                            
                                exp1 = exp(1i*k*dot(d1,y-x));
                                exp2 = exp(1i*k*dot(d2,y-x));
                                Psi2 = 1i*(exp1-exp2)/k;
                                Psi1 = (exp1+exp2)/2;
                                dPsi2dny = dot(d2,ny)*exp2 - dot(d1,ny)*exp1;
                                dPsi1dny = 1i*k*(dot(d1,ny)*exp1+dot(d2,ny)*exp2)/2;

                                dPhi_kTemp = dPhi_kdny(xmy,r,ny);
                                Phi_kTemp = Phi_k(r);

                                Psi1_integral     = Psi1_integral    + Psi1*dPhi_kTemp*fact_y; 
                                Psi2_integral     = Psi2_integral    + Psi2*dPhi_kTemp*fact_y; 
                                Psi3_integral     = Psi3_integral    + (dPsi3dny*Phi_kTemp-Psi3*dPhi_kTemp)*fact_y; 
                                Psi4_integral     = Psi4_integral    + (dPsi4dny*Phi_kTemp-Psi4*dPhi_kTemp)*fact_y; 
                                dPsi1dny_integral = dPsi1dny_integral + dPsi1dny*Phi_kTemp*fact_y;
                                dPsi2dny_integral = dPsi2dny_integral + dPsi2dny*Phi_kTemp*fact_y;

                                CBIE = CBIE + dPhi_kTemp*R_y*fact_y;
                            end
                        end
                    end
                end
            end
            idxCol(:,e_y) = sctr_y;
            if useCBIE
                A_e_temp(:,:,e_y) = A_e_temp(:,:,e_y) + R_x.'*CBIE*fact_x;
            end
        end
        if useCBIE            
            A_e_temp(:,:,e_x) = A_e_temp(:,:,e_x) + R_x.'*R_x*(dPsi1dny_integral - Psi1_integral - 1)*fact_x;
            F_e = F_e - R_x.'*p_inc(x).'*fact_x;
        end
    end
    
    idxRow(:,e_x) = sctr_x.';
    Avalues(:,:,e_x) = matrixAssembly(A_e_temp, idxCol, n_en, noDofs, noElems, 1);
    Fvalues(:,e_x,:) = F_e;
end
A = matrixAssembly(Avalues, idxRow, n_en, noDofs, noElems, 2);
FF = zeros(noDofs,no_angles);
parfor i = 1:no_angles
    FF(:,i) = vectorAssembly(Fvalues(:,:,i), idxRow, noDofs);
end


