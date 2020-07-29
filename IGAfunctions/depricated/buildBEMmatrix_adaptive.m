function [A, FF] = buildBEMmatrix_adaptive(varCol)
error('Depricated. Use buildCBEMmatrix instead')

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;

noElems = varCol.noElems;
index = varCol.index;
nurbs = varCol.nurbs;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
n_xi = varCol.nurbs.number(1);
n_eta = varCol.nurbs.number(2);

p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

uniqueXi = unique(Xi);
uniqueEta = unique(Eta);
noElementsXi = length(uniqueXi)-1;
noElementsEta = length(uniqueEta)-1;


dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

model = varCol.model;
k = varCol.k;
alpha = 1i/k;

Phi_0 = @(r)           1/(4*pi*r);
Phi_k = @(r) exp(1i*k*r)/(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)/r^2*             (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)/r^2*(1 - 1i*k*r)*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

h_min = findMaxElementDiameter(noElems, nurbs);

switch varCol.formulation(2:end)
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



if strcmp(model, 'PS')
    no_angles = 1;
    P_inc = NaN;
    dpdn = varCol.dpdn;
    dP_inc = NaN;
    analytic = varCol.analytic;
    gAnalytic = varCol.gAnalytic;
else
    no_angles = length(varCol.alpha_s);
    P_inc = varCol.P_inc;
    dP_inc = varCol.dP_inc;
    analytic = @(v) varCol.analytic(v) + P_inc(v);
    gAnalytic = @(v) varCol.gAnalytic(v) + varCol.gP_inc(v);
    dpdn = 0;
end
exteriorProblem = true;
if exteriorProblem
    sgn = -1;
else
    sgn = 1;
end

%% Create collocation points based on the Greville abscissae
if false
    [cg_xi, grev_xi] = CauchyGalerkin(p_xi, n_xi, Xi);
    [cg_eta, grev_eta] = CauchyGalerkin(p_eta, n_eta, Eta);

    n_cp = noDofs - length(dofsToRemove);

    cp_p = zeros(n_cp,2);
    counter = 1;
    counter2 = 1;
    for j = 1:n_eta
        eta_bar = cg_eta(j);
        for i = 1:n_xi
            if ~any(dofsToRemove == counter)
                xi_bar = cg_xi(i);      
                cp_p(counter2,:) = [xi_bar, eta_bar]; 
                counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
else
    n_cp = noDofs - length(dofsToRemove);
    switch varCol.model
        case {'SS', 'S1', 'S2', 'S3', 'PS', 'MS', 'EL'}
            eps_greville_xi = 0;
            eps_greville_eta = 0;
%             eps_greville_xi = 1e-1*min(uniqueXi(2:end)-uniqueXi(1:end-1));
%             eps_greville_eta = 1e-1*min(uniqueEta(2:end)-uniqueEta(1:end-1));
        otherwise
            eps_greville_xi = 1e-1*min(uniqueXi(2:end)-uniqueXi(1:end-1));
            eps_greville_eta = 1e-1*min(uniqueEta(2:end)-uniqueEta(1:end-1));
    end
    cp_p = zeros(n_cp,2);
    cp = zeros(n_cp,3);
    counter = 1;
    counter2 = 1;
    for j = 1:n_eta
        eta_bar = sum(Eta(j+1:j+p_eta))/p_eta;
        for i = 1:n_xi
            if ~any(dofsToRemove == counter)
                xi_bar = sum(Xi(i+1:i+p_xi))/p_xi;
                if ismember(xi_bar, Xi)
                    if Xi(i+1) == Xi(i+p_xi+1)
                        xi_bar = xi_bar - eps_greville_xi;
                    else
                        xi_bar = xi_bar + eps_greville_xi;
                    end
                end
                if ismember(eta_bar, Eta)
                    if Eta(j+1) == Eta(j+p_eta+1)
                        eta_bar = eta_bar - eps_greville_eta;
                    else
                        eta_bar = eta_bar + eps_greville_eta;
                    end
                end

                cp_p(counter2,:) = [xi_bar, eta_bar]; 
                cp(counter2,:) = evaluateNURBS(nurbs, [xi_bar, eta_bar]);
                counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
end


%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));

errorAvg = 0;

computeLHS = 1;
plotMesh = 1;
for i = 520
    
% parfor i = 1:n_cp
    totArea = 0;
    Phi_k_integralExp = 0;
    d2Phi_kdnxdny_integrals = 0;
    dPhi_0dny_integral = 0;
    d2Phi_0dnxdny_integral = 0;
    ugly_integral = zeros(3,1);
    A_row = zeros(1, noDofs);
    lhs = 0;
    
    xi_x = cp_p(i,1);
    eta_x = cp_p(i,2);
%     if xi_x > 0.2 && xi_x < 0.8 && eta_x > 0.2 && eta_x < 0.8
%         plotMesh = 1;
%     else
%         plotMesh = 0;
%     end
    xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
    eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
    e_x = xi_idx + noElementsXi*(eta_idx-1);
    sctr_x = element(e_x,:);
%     hold on
%     plot3(x(1),x(2),x(3), '*')
%     drawnow
    if useHBIE
        [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, weights);
        pts_x = controlPts(sctr_x,:);
        x = R_x*pts_x;
        J_temp = [dR_xdxi; dR_xdeta]*pts_x;
        m_1 = J_temp(1,:);
        m_2 = J_temp(2,:);
        crossProd_x = cross(m_1,m_2);
        h_xi = norm(m_1);
        h_eta = norm(m_2);
        
        if (eta_x == 0 || eta_x == 1) && (strcmp(model,'SS') || strcmp(model,'S1') || strcmp(model,'S2') ...
                                            || strcmp(model,'PS') || strcmp(model,'MS') || strcmp(model,'EL'))
            v_2 = m_2/h_eta;
            nx = x.'/norm(x);
            v_3 = nx.';
            v_1 = cross(v_2,v_3);
            J_x = [1, 0; 0, 1/h_eta];
        else
            e_xi = m_1/h_xi;
            e_eta = m_2/h_eta;
            v_1 = m_1/norm(m_1);
            nx = crossProd_x'/norm(crossProd_x);
            v_3 = nx.';
            v_2 = cross(v_3,v_1);
            cosT = dot(e_xi,e_eta);
            sinT = dot(v_2,e_eta);
            J_x = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
        end
    else
        R_x = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, weights);
        pts_x = controlPts(sctr_x,:);
        x = R_x*pts_x;
        nx = NaN;        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotMesh
        close all
        hold off
%         plotNURBS(varCol.nurbs,[40 40], 1, getColor(1), 1);

        varCol.x = x;
        varCol.colorFun = @(xmy,r,nx,ny) Phi_k(r);
%         varCol_temp.colorFun = @(xmy,r,nx,ny) 1/r^2;
%         varCol_temp.colorFun = @(xmy,r,nx,ny) d2Phi_kdnxdny(xmy,r,nx,ny);
        varCol.xi_x = xi_x;
        varCol.eta_x = eta_x;
        varCol.h_min = h_min;
        if true
            plotNURBS(nurbs, [40,40], 1, NaN, NaN, NaN, varCol);
        else
            plotNURBSadaptive(nurbs, [40,40], 1, NaN, varCol);
%         camlight
    %         colormap jet
    %         savefig(['fig_i_' num2str(i) '_woMesh.fig'])  
        end
        axis equal
        axis off
        set(gca, 'Color', 'none');
        view(78,22)
%         colorbar
        drawnow
        hold on
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FF_temp = zeros(1, no_angles);

    for e = 1:noElems        
        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);


        sctr = element(e,:);
%         sctr2 = element2(e,:);
%         R_x_e = R_x(sctr2);

        pts = controlPts(sctr,:);
        CBIE = zeros(1, n_en);
        HBIE = zeros(1, n_en);
        

        [FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, d2Phi_kdnxdny_integrals] ...
            = adaptiveIntegration(FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, ...
            d2Phi_kdnxdny_integrals, Xi_e, Eta_e, p_xi, p_eta, Xi, Eta, weights, pts, nurbs, computeLHS, k, useCBIE, useHBIE, model, x, ...
            xi_x, eta_x, h_min, analytic, gAnalytic, alpha, nx, plotMesh, dpdn);
        if useCBIE
            for j = 1:n_en
                A_row(sctr(j)) = A_row(sctr(j)) + CBIE(j);
            end
        end
        if useHBIE
            for j = 1:n_en
                A_row(sctr(j)) = A_row(sctr(j)) + HBIE(j);
            end
        end
    end
    if useCBIE
        for j = 1:n_en
            A_row(sctr_x(j)) = A_row(sctr_x(j)) - R_x(j)*(0.5*(1-sgn) + dPhi_0dny_integral);
        end
        if computeLHS
            lhs = lhs - analytic(x)*(0.5*(1-sgn) + dPhi_0dny_integral);
        end
    end
    if useHBIE
        dphidv = J_x*[dR_xdxi; dR_xdeta];
%         temp = (J_x(1,:)*(v_1*ugly_integral) + J_x(2,:)*(v_2*ugly_integral))*[dR_xdxi; dR_xdeta];
        gR_x = dphidv(1,:).'*v_1 + dphidv(2,:).'*v_2;
        temp = gR_x*ugly_integral;
        for j = 1:n_en
            A_row(sctr_x(j)) = A_row(sctr_x(j)) + alpha*(temp(j) - R_x(j)*d2Phi_0dnxdny_integral);
        end
        if computeLHS
            lhs = lhs + alpha*(gAnalytic(x)*ugly_integral - analytic(x)*d2Phi_0dnxdny_integral);
        end
    end
    A(i,:) = A_row;
    if strcmp(model,'PS')
        FF(i,:) = FF_temp;
        if useHBIE
            FF(i,:) = FF(i,:) + alpha*dpdn(x,nx)*(dPhi_0dny_integral + 0.5*(1-sgn) - v_3*ugly_integral);
%             FF(i,:) = FF(i,:) + alpha*dpdn(x,nx)*(dPhi_0dny_integral + 0.5*(1-sgn));
        end
    else
        if useCBIE
            FF(i,:) = -P_inc(x).';
        end
        if useHBIE
            FF(i,:) = FF(i,:) - alpha*dP_inc(x,nx).';
        end
    end
    if computeLHS
        if FF(i,:) ~= 0
            errorAvg = errorAvg + abs(lhs-FF(i,:))/max(1,abs(FF(i,:)))
        else
            errorAvg = errorAvg + abs(lhs)
        end
    end
    keyboard
    if plotMesh
%         savefig(['fig_i_' num2str(i) '.fig'])
%         keyboard
    end
%     keyboard
%     keyboard
%     if useHBIE
%         abs(Phi_k_integralExp-d2Phi_kdnxdny_integrals)/abs(d2Phi_kdnxdny_integrals)
%     end
%     plot3(x(1),x(2),x(3), '*')
%     keyboard
%     dPhi_0dny_integral+0.5
%     d2Phi_0dnxdny_integral
%     if strcmp(varCol.model,'SS') || strcmp(varCol.model,'PS')
%         R_o = 0.5;
%         errorInTotArea = abs((totArea-4*pi*R_o^2)/(4*pi*R_o^2))
% %         keyboard
%     elseif strcmp(varCol.model,'PL')
%         errorInTotArea = abs(totArea-1)
%     end
end

if computeLHS
    100*errorAvg/n_cp
end

function [FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, d2Phi_kdnxdny_integrals] ...
    = adaptiveIntegration(FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, d2Phi_kdnxdny_integrals, ...
    Xi_e, Eta_e, p_xi, p_eta, Xi, Eta, weights, pts, nurbs, computeLHS, k, useCBIE, useHBIE, model, x, ...
    xi_x, eta_x, h_min, analytic, gAnalytic, alpha, nx, plotMesh, dpdn)
runSimpsonAdaptivity = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotMesh
    hold on
    plot3(x(1),x(2),x(3),'*','color', 'red')
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_0 = @(r)           1/(4*pi*r);
Phi_k = @(r) exp(1i*k*r)/(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)/r^2*             (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)/r^2*(1 - 1i*k*r)*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

x_1 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(1)+eps]).';
x_2 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(1)+eps]).';
x_3 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(2)-eps]).';
x_4 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(2)-eps]).';


l = min([norm(x-x_1), norm(x-x_2), norm(x-x_3), norm(x-x_4)]);
h_1 = norm(x_1-x_3);
h_2 = norm(x_2-x_4);
h = max(h_1,h_2);
p_max = max(p_xi,p_eta);
coeff2 = 0.3; % scale the inverse density of GPs
coeff3 = 10; % measure of the size of element containing singularity
coeff = 1; %mesh 1 -> 5 er optimal
coeff4 = 3;% round((coeff2/coeff3)^2);
n_gp_xi = coeff*(p_xi+1);
n_gp_eta = coeff*(p_eta+1);
n_gp_rho = coeff4*(p_max+1);
n_gp_theta = coeff4*(p_max+1);
xi_h = max(norm(x_1-x_2),norm(x_3-x_4));
eta_h = max(norm(x_2-x_3),norm(x_1-x_4));

if eta_x == 0 || eta_x == 1
    xi_x = Xi_e(1); 
elseif xi_x == 0 && Xi_e(2) == 1
    xi_x = 1;
elseif xi_x == 1 && Xi_e(1) == 0 % probably redundant
    xi_x = 0;
end

singularityIsInElement = (Xi_e(1) <= xi_x+10*eps) && (xi_x-10*eps <= Xi_e(2)) && (Eta_e(1) <= eta_x+10*eps) && (eta_x-10*eps <= Eta_e(2));
% singularityIsInElement = Xi_e(1) <= xi_x && xi_x <= Xi_e(2) && Eta_e(1) <= eta_x && eta_x <= Eta_e(2);
if ~runSimpsonAdaptivity
    if (singularityIsInElement && h/h_min > coeff3) || (~singularityIsInElement && h/l > coeff2)
%         keyboard
        if xi_h/eta_h > 2
            n_div_xi = round(xi_h/eta_h);
            n_div_eta = 1;
        elseif eta_h/xi_h > 2
            n_div_xi = 1;
            n_div_eta = round(eta_h/xi_h);
        elseif eta_h > xi_h
            n_div_xi = 2;
            n_div_eta = round(2*eta_h/xi_h);
        else
            n_div_xi = round(2*xi_h/eta_h);
            n_div_eta = 2;
        end

        Xi_e_arr  = linspace(Xi_e(1),Xi_e(2),n_div_xi+1);
        Eta_e_arr = linspace(Eta_e(1),Eta_e(2),n_div_eta+1);
        if singularityIsInElement
            [spn, I] = min(abs(Xi_e_arr(2:end-1) - xi_x));
            I = I + 1;
            if spn < 1/12*(Xi_e(2)-Xi_e(1))/n_div_xi
                Xi_e_arr = [Xi_e_arr(1:I-1), 1/3*Xi_e_arr(I+1)+2/3*Xi_e_arr(I-1), 2/3*Xi_e_arr(I+1)+1/3*Xi_e_arr(I-1), Xi_e_arr(I+1:end)];
                n_div_xi = n_div_xi + 1;
            end
            [spn, I] = min(abs(Eta_e_arr(2:end-1) - eta_x));
            I = I + 1;
            if spn < 1/12*(Eta_e(2)-Eta_e(1))/n_div_eta
                Eta_e_arr = [Eta_e_arr(1:I-1), 1/3*Eta_e_arr(I+1)+2/3*Eta_e_arr(I-1), 2/3*Eta_e_arr(I+1)+1/3*Eta_e_arr(I-1), Eta_e_arr(I+1:end)];
                n_div_eta = n_div_eta + 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plotMesh
            nn = 100;
            EEta = linspace(Eta_e(1),Eta_e(2),nn);
            for ii = 1:length(Xi_e_arr)
                if ismember(Xi_e_arr(ii),Xi)
                    continue
                end
                XX = zeros(nn,3);
                for jj = 1:nn
                    XX(jj,:) = evaluateNURBS(nurbs, [Xi_e_arr(ii), EEta(jj)]);
                end
                hold on
                plot3(XX(:,1),XX(:,2),XX(:,3),'white')
                hold off
            end
            XXi = linspace(Xi_e(1),Xi_e(2),nn);
            for ii = 1:length(Eta_e_arr)
                if ismember(Eta_e_arr(ii),Eta)
                    continue
                end
                XX = zeros(nn,3);
                for jj = 1:nn
                    XX(jj,:) = evaluateNURBS(nurbs, [XXi(jj), Eta_e_arr(ii)]);
                end
                hold on
                plot3(XX(:,1),XX(:,2),XX(:,3),'white')
                hold off
            end
            drawnow
        end
%         keyboard
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i_eta = 1:n_div_eta
            Eta_e_sub = Eta_e_arr(i_eta:i_eta+1);
            for i_xi = 1:n_div_xi
                Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);
                [FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, d2Phi_kdnxdny_integrals] ...
                    = adaptiveIntegration(FF_temp, dPhi_0dny_integral, CBIE, lhs, d2Phi_0dnxdny_integral, HBIE, ugly_integral, totArea, Phi_k_integralExp, ...
                    d2Phi_kdnxdny_integrals, Xi_e_sub, Eta_e_sub, p_xi, p_eta, Xi, Eta, weights, pts, nurbs, computeLHS, k, useCBIE, useHBIE, model, x, ...
                    xi_x, eta_x, h_min, analytic, gAnalytic, alpha, nx, plotMesh, dpdn);
            end
        end
    else      
        if singularityIsInElement
            xi_x_t = parametric2parentSpace(Xi_e, xi_x);
            eta_x_t = parametric2parentSpace(Eta_e, eta_x);
            theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
            theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
            theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
            theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            for area = {'South', 'East', 'North', 'West'}
                switch area{1}
                    case 'South'
                        if abs(eta_x - Eta_e(1)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x3 theta_x4];
                    case 'East'
                        if abs(xi_x - Xi_e(2)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x4 theta_x1];
                    case 'North'
                        if abs(eta_x - Eta_e(2)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x1 theta_x2];
                    case 'West'
                        if abs(xi_x - Xi_e(1)) < 10*eps
                            continue
                        end
                        if theta_x3 < 0
                            thetaRange = [theta_x2 theta_x3+2*pi];
                        else
                            thetaRange = [theta_x2 theta_x3];
                        end
                end

                [W2D,Q2D] = gaussianQuadNURBS(n_gp_rho,n_gp_theta);

                for gp = 1:size(W2D,1)
                    pt = Q2D(gp,:);
                    wt = W2D(gp);

                    rho_t = parent2ParametricSpace([0, 1],   pt(1));
                    theta = parent2ParametricSpace(thetaRange,pt(2));
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

                    xi_t  = xi_x_t + rho*cos(theta);
                    eta_t = eta_x_t + rho*sin(theta);

                    J_3 = rho;
                    J_4 = rho_hat;
                    J_5 = 0.25*(thetaRange(2)-thetaRange(1));

                    xi = parent2ParametricSpace(Xi_e, xi_t);
                    eta = parent2ParametricSpace(Eta_e, eta_t);

                    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                    J = [dR_ydxi; dR_ydeta]*pts;
                    crossProd = cross(J(1,:),J(2,:));
                    J_1 = norm(crossProd);
                    ny = crossProd'/J_1;

                    y = R_y*pts;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 hold on
    %                 plot3(y(1),y(2),y(3),'*')
    %                 hold off
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xmy = x-y;
                    r = norm(xmy);
                    fact = J_1*J_2*J_3*J_4*J_5*wt;


                    if strcmp(model,'PS')
                        if useCBIE
                            FF_temp = FF_temp + Phi_k(r)*dpdn(y,ny)*fact;
                        end
                        if useHBIE
                            FF_temp = FF_temp + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                        end
                    end
                    dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                    if useCBIE
                        CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                        if computeLHS
                            lhs = lhs + dPhi_kdny(xmy,r,ny)*analytic(y)*fact;
                        end
                    end

                    if useHBIE
                        d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                        HBIE = HBIE + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                        ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                            + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                        if computeLHS
                            lhs = lhs + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*analytic(y)*fact;
                        end
                    end
                    totArea = totArea + fact; 
                    if useHBIE
                        Phi_k_integralExp = Phi_k_integralExp + k^2*Phi_k(r)*dot(nx,ny)*fact;
                        d2Phi_kdnxdny_integrals = d2Phi_kdnxdny_integrals + (d2Phi_kdnxdny(xmy,r,nx,ny) - d2Phi_0dnxdny(xmy,r,nx,ny))*fact;
                    end
                end  
            end
        else
            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
            [W2D,Q2D] = gaussianQuadNURBS(n_gp_xi,n_gp_eta);
            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);
                wt = W2D(gp);

                xi  = parent2ParametricSpace(Xi_e, pt(1));
                eta = parent2ParametricSpace(Eta_e,pt(2));
                [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                J = [dR_ydxi; dR_ydeta]*pts;
                crossProd = cross(J(1,:),J(2,:));
                J_1 = norm(crossProd);
                ny = crossProd'/J_1;

                y = R_y*pts;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             hold on
    %             plot3(y(1),y(2),y(3),'*')
    %             hold off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                xmy = x-y;
                r = norm(xmy);
                fact = J_1*J_2*wt;

                if strcmp(model,'PS')
                    if useCBIE
                        FF_temp = FF_temp + Phi_k(r)*dpdn(y,ny)*fact;
                    end
                    if useHBIE
                        FF_temp = FF_temp + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                    end
                end
                dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                if useCBIE
                    CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                    if computeLHS
                        lhs = lhs + dPhi_kdny(xmy,r,ny)*analytic(y)*fact;
                    end
                end
                if useHBIE
                    d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                    HBIE = HBIE + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                    ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                        + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                    if computeLHS
                        lhs = lhs + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*analytic(y)*fact;
                    end
                end
                totArea = totArea + fact;
                if useHBIE
                    Phi_k_integralExp = Phi_k_integralExp + k^2*Phi_k(r)*dot(nx,ny)*fact;
                    d2Phi_kdnxdny_integrals = d2Phi_kdnxdny_integrals + (d2Phi_kdnxdny(xmy,r,nx,ny) - d2Phi_0dnxdny(xmy,r,nx,ny))*fact;
                end
            end
        end
    end
else
%     keyboard
    if runSimpsonAdaptivity
        if singularityIsInElement
            xi_x_t = parametric2parentSpace(Xi_e, xi_x);
            eta_x_t = parametric2parentSpace(Eta_e, eta_x);
            theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
            theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
            theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
            theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            for area = {'South', 'East', 'North', 'West'}
                switch area{1}
                    case 'South'
                        if abs(eta_x - Eta_e(1)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x3 theta_x4];
                    case 'East'
                        if abs(xi_x - Xi_e(2)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x4 theta_x1];
                    case 'North'
                        if abs(eta_x - Eta_e(2)) < 10*eps
                            continue
                        end
                        thetaRange = [theta_x1 theta_x2];
                    case 'West'
                        if abs(xi_x - Xi_e(1)) < 10*eps
                            continue
                        end
                        if theta_x3 < 0
                            thetaRange = [theta_x2 theta_x3+2*pi];
                        else
                            thetaRange = [theta_x2 theta_x3];
                        end
                end

                [W2D,Q2D] = gaussianQuadNURBS(n_gp_rho,n_gp_theta);

                for gp = 1:size(W2D,1)
                    pt = Q2D(gp,:);
                    wt = W2D(gp);

                    rho_t = parent2ParametricSpace([0, 1],   pt(1));
                    theta = parent2ParametricSpace(thetaRange,pt(2));
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

                    xi_t  = xi_x_t + rho*cos(theta);
                    eta_t = eta_x_t + rho*sin(theta);

                    J_3 = rho;
                    J_4 = rho_hat;
                    J_5 = 0.25*(thetaRange(2)-thetaRange(1));

                    xi = parent2ParametricSpace(Xi_e, xi_t);
                    eta = parent2ParametricSpace(Eta_e, eta_t);

                    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                    J = [dR_ydxi; dR_ydeta]*pts;
                    crossProd = cross(J(1,:),J(2,:));
                    J_1 = norm(crossProd);
                    ny = crossProd'/J_1;

                    y = R_y*pts;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 hold on
    %                 plot3(y(1),y(2),y(3),'*')
    %                 hold off
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xmy = x-y;
                    r = norm(xmy);
                    fact = J_1*J_2*J_3*J_4*J_5*wt;


                    if strcmp(model,'PS')
                        if useCBIE
                            FF_temp = FF_temp + Phi_k(r)*dpdn(y,ny)*fact;
                        end
                        if useHBIE
                            FF_temp = FF_temp + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                        end
                    end
                    dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                    if useCBIE
                        CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                        if computeLHS
                            lhs = lhs + dPhi_kdny(xmy,r,ny)*analytic(y)*fact;
                        end
                    end

                    if useHBIE
                        d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                        HBIE = HBIE + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                        ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                            + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                        if computeLHS
                            lhs = lhs + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*analytic(y)*fact;
                        end
                    end
                    totArea = totArea + fact; 
                    if useHBIE
                        Phi_k_integralExp = Phi_k_integralExp + k^2*Phi_k(r)*dot(nx,ny)*fact;
                        d2Phi_kdnxdny_integrals = d2Phi_kdnxdny_integrals + (d2Phi_kdnxdny(xmy,r,nx,ny) - d2Phi_0dnxdny(xmy,r,nx,ny))*fact;
                    end
                end  
            end
        else
            x_5 = evaluateNURBS(nurbs, [mean(Xi_e),  mean(Eta_e)]).';
            
            l = norm(x-x_5);
            h_1 = norm(x_1-x_3);
            h_2 = norm(x_2-x_4);
            h = max(h_1,h_2);
            n_div = round(2*h/l + 1);
%             n_div = 1;
            Xi_e_arr  = linspace(Xi_e(1),Xi_e(2),n_div+1);
            Eta_e_arr = linspace(Eta_e(1),Eta_e(2),n_div+1);
            for i_eta = 1:n_div
                Eta_e_sub = Eta_e_arr(i_eta:i_eta+1);
                for i_xi = 1:n_div
                    Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);
                    J_2 = 0.25*(Xi_e_sub(2)-Xi_e_sub(1))*(Eta_e_sub(2)-Eta_e_sub(1));
                    if plotMesh
                        nn = 100;
                        EEta = linspace(Eta_e_sub(1),Eta_e_sub(2),nn);
                        for ii = 1:length(Xi_e_arr)
                            if ismember(Xi_e_arr(ii),Xi)
                                continue
                            end
                            XX = zeros(nn,3);
                            for jj = 1:nn
                                XX(jj,:) = evaluateNURBS(nurbs, [Xi_e_arr(ii), EEta(jj)]);
                            end
                            hold on
                            plot3(XX(:,1),XX(:,2),XX(:,3),'white')
                            hold off
                        end
                        XXi = linspace(Xi_e_sub(1),Xi_e_sub(2),nn);
                        for ii = 1:length(Eta_e_arr)
                            if ismember(Eta_e_arr(ii),Eta)
                                continue
                            end
                            XX = zeros(nn,3);
                            for jj = 1:nn
                                XX(jj,:) = evaluateNURBS(nurbs, [XXi(jj), Eta_e_arr(ii)]);
                            end
                            hold on
                            plot3(XX(:,1),XX(:,2),XX(:,3),'white')
                            hold off
                        end
                        drawnow
                    end
                    [W2D,Q2D] = gaussianQuadNURBS(n_gp_xi,n_gp_eta);
                    for gp = 1:size(W2D,1)
                        pt = Q2D(gp,:);
                        wt = W2D(gp);

                        xi  = parent2ParametricSpace(Xi_e_sub, pt(1));
                        eta = parent2ParametricSpace(Eta_e_sub,pt(2));
                        [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                        J = [dR_ydxi; dR_ydeta]*pts;
                        crossProd = cross(J(1,:),J(2,:));
                        J_1 = norm(crossProd);
                        ny = crossProd'/J_1;

                        y = R_y*pts;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             hold on
            %             plot3(y(1),y(2),y(3),'*')
            %             hold off
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        xmy = x-y;
                        r = norm(xmy);
                        fact = J_1*J_2*wt;

                        if strcmp(model,'PS')
                            if useCBIE
                                FF_temp = FF_temp + Phi_k(r)*dpdn(y,ny)*fact;
                            end
                            if useHBIE
                                FF_temp = FF_temp + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                            end
                        end
                        dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                        if useCBIE
                            CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                            if computeLHS
                                lhs = lhs + dPhi_kdny(xmy,r,ny)*analytic(y)*fact;
                            end
                        end
                        if useHBIE
                            d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                            HBIE = HBIE + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                            ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                                + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                            if computeLHS
                                lhs = lhs + alpha*d2Phi_kdnxdny(xmy,r,nx,ny)*analytic(y)*fact;
                            end
                        end
                        totArea = totArea + fact;
                        if useHBIE
                            Phi_k_integralExp = Phi_k_integralExp + k^2*Phi_k(r)*dot(nx,ny)*fact;
                            d2Phi_kdnxdny_integrals = d2Phi_kdnxdny_integrals + (d2Phi_kdnxdny(xmy,r,nx,ny) - d2Phi_0dnxdny(xmy,r,nx,ny))*fact;
                        end
                    end
                end
            end
        end
    end
end