function R_h = evaluateBEMerrorFunc(varCol, xi_x, eta_x, U)

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


Phi_0 = @(r)           1/(4*pi*r);
Phi_k = @(r) exp(1i*k*r)/(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)/r^2*             (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)/r^2*(1 - 1i*k*r)*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

if strcmp(model, 'PS')
    no_angles = 1;
    P_inc = NaN;
    dpdn = varCol.dpdn;
    dP_inc = NaN;
else
    no_angles = length(varCol.alpha_s);
    P_inc = varCol.P_inc;
    dP_inc = varCol.dP_inc;
    dpdn = 0;
end
exteriorProblem = true;
if exteriorProblem
    sgn = -1;
else
    sgn = 1;
end



%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1);
p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(3*p_max+1,3*p_max+1);

dPhi_0dny_integral = 0;
d2Phi_0dnxdny_integral = 0;
ugly_integral = zeros(3,1);
A_row = zeros(1, noDofs);

xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
e_x = xi_idx + noElementsXi*(eta_idx-1);
sctr_x = element(e_x,:);

if useHBIE
    [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, weights);
    pts_x = controlPts(sctr_x,:);
    J_temp = [dR_xdxi; dR_xdeta]*pts_x;
    m_1 = J_temp(1,:);
    m_2 = J_temp(2,:);
    crossProd_x = cross(m_1,m_2);
    h_xi = norm(m_1);
    h_eta = norm(m_2);
    e_xi = m_1/h_xi;
    e_eta = m_2/h_eta;
    x = R_x*controlPts(sctr_x,:);

    if (eta_x == 0 || eta_x == 1) && (strcmp(model,'SS') || strcmp(model,'S1') || strcmp(model,'S2') ...
                                        || strcmp(model,'PS') || strcmp(model,'MS') || strcmp(model,'EL'))
        v_2 = m_2/h_eta;
        nx = x.'/norm(x);
        v_3 = nx.';
        v_1 = cross(v_2,v_3);
        J_x = [1, 0; 0, 1/h_eta];
    else
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
    x = R_x*controlPts(sctr_x,:);
    nx = NaN;        
end


F = 0;

for e = 1:noElems        
    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);


    sctr = element(e,:);

    pts = controlPts(sctr,:);
    if useCBIE
        CBIE = zeros(1, n_en);
    end
    if useHBIE
        HBIE = zeros(1, n_en);
    end
    if eta_x == 0
        xi_x = Xi_e(1);
    elseif eta_x == 1 || (xi_x == 0 && Xi_e(2) == 1)
        xi_x = Xi_e(2);
    elseif xi_x == 1 && Xi_e(1) == 0
        xi_x = 0;
    end
    if Xi_e(1) <= xi_x && xi_x <= Xi_e(2) && Eta_e(1) <= eta_x && eta_x <= Eta_e(2)
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
            for gp = 1:size(W2D_2,1)
                pt = Q2D_2(gp,:);
                wt = W2D_2(gp);

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

                xi = parent2ParametricSpace(Xi_e, xi_t);
                eta = parent2ParametricSpace(Eta_e, eta_t);

                [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                J = pts'*[dR_ydxi' dR_ydeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                ny = crossProd/J_1;

                J_3 = rho;
                J_4 = rho_hat;
                J_5 = 0.25*(thetaRange(2)-thetaRange(1));

                y = R_y*pts;
                
                xmy = x-y;
                r = norm(xmy);
                fact = J_1*J_2*J_3*J_4*J_5*wt;


                if strcmp(model,'PS')
                    if useCBIE
                        F = F + Phi_k(r)*dpdn(y,ny)*fact;
                    end
                    if useHBIE
                        F = F + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                    end
                end
                dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                if useCBIE
                    CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                end
                if useHBIE
                    d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                    HBIE = HBIE + d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                    ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                        + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                end
%                     totArea = totArea + fact;    
            end  
        end
    else
        x_1 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(1)+eps]).';
        x_2 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(1)+eps]).';
        x_3 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(2)-eps]).';
        x_4 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(2)-eps]).';
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
                for gp = 1:size(W2D,1)
                    pt = Q2D(gp,:);
                    wt = W2D(gp);

                    xi  = parent2ParametricSpace(Xi_e_sub, pt(1));
                    eta = parent2ParametricSpace(Eta_e_sub,pt(2));
                    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);

                    J = pts'*[dR_ydxi' dR_ydeta'];
                    crossProd = cross(J(:,1),J(:,2));
                    J_1 = norm(crossProd);
                    ny = crossProd/J_1;

                    y = R_y*pts;
                    
                    xmy = x-y;
                    r = norm(xmy);
                    fact = J_1*J_2*wt;

                    if strcmp(model,'PS')
                        if useCBIE
                            F = F + Phi_k(r)*dpdn(y,ny)*fact;
                        end
                        if useHBIE
                            F = F + alpha*dPhi_kdnx(xmy,r,nx)*dpdn(y,ny)*fact;
                        end
                    end
                    dPhi_0dny_integral = dPhi_0dny_integral + dPhi_0dny(xmy,r,ny)*fact; 
                    if useCBIE
                        CBIE = CBIE + dPhi_kdny(xmy,r,ny)*R_y*fact;
                    end
                    if useHBIE
                        d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + d2Phi_0dnxdny(xmy,r,nx,ny)*fact;
                        HBIE = HBIE + d2Phi_kdnxdny(xmy,r,nx,ny)*R_y*fact;
                        ugly_integral = ugly_integral + (dPhi_0dnx(xmy,r,nx)*ny + dPhi_0dny(xmy,r,ny)*nx ...
                                                                            + d2Phi_0dnxdny(xmy,r,nx,ny)*xmy.')*fact;
                    end
%                         totArea = totArea + fact;
                end
            end
        end
    end
    if useCBIE
        for j = 1:n_en
            A_row(sctr(j)) = A_row(sctr(j)) + CBIE(j);
        end
    end
    if useHBIE
        for j = 1:n_en
            A_row(sctr(j)) = A_row(sctr(j)) + alpha*HBIE(j);
        end
    end
end

if useCBIE
    for j = 1:n_en
        A_row(sctr_x(j)) = A_row(sctr_x(j)) - R_x(j)*(0.5*(1-sgn) + dPhi_0dny_integral);
    end
end
if useHBIE
    temp = (J_x(1,:)*(v_1*ugly_integral) + J_x(2,:)*(v_2*ugly_integral))*[dR_xdxi; dR_xdeta];
    for j = 1:n_en
        A_row(sctr_x(j)) = A_row(sctr_x(j)) + alpha*(-R_x(j)*d2Phi_0dnxdny_integral + temp(j));
    end
end
if strcmp(model,'PS')
    if useHBIE
        F = F + alpha*dpdn(x,nx)*(dPhi_0dny_integral + 0.5*(1-sgn) - v_3*ugly_integral);
    end
else
    if useCBIE
        F = F - P_inc(x).';
    end
    if useHBIE
        F = F - alpha*dP_inc(x,nx).';
    end
end
R_h = A_row*U - F;

