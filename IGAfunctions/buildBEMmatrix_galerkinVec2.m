function [A, FF, varCol] = buildBEMmatrix_galerkinVec2(varCol)

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
patches = varCol.patches;

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
noDofs = varCol.noDofs;
agpBEM = varCol.agpBEM;
exteriorProblem = varCol.exteriorProblem;

quadMethodBEMsimpson = strcmp(varCol.quadMethodBEM,'Simpson');

Eps = 10*eps;

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

Phi_0 = @(r)           1./(4*pi*r);
Phi_k = @(r) exp(1i*k*r)./(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)./r.^2.*             sum(xmy.*ny,2);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)./r.^2.*(1 - 1i*k*r).*sum(xmy.*ny,2);

dPhi_0dnx = @(xmy,r,nx) -Phi_0(r)./r.^2.*             (xmy*nx);
dPhi_kdnx = @(xmy,r,nx) -Phi_k(r)./r.^2.*(1 - 1i*k*r).*(xmy*nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)./r.^2.*((ny*nx)           - 3./r.^2                .*(xmy*nx).*sum(xmy.*ny,2));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)./r.^2.*((ny*nx).*(1-1i*k*r)+(k^2+3./r.^2.*(1i*k*r-1)).*(xmy*nx).*sum(xmy.*ny,2));

SHBC = strcmp(varCol.BC, 'SHBC');
if SHBC
    no_angles = length(varCol.alpha_s);
    p_inc = varCol.p_inc;
    dp_inc = varCol.dp_inc;
else
    no_angles = 1;
    p_inc = NaN;
    dp_inc = NaN;
end
dpdn = varCol.dpdn;

if exteriorProblem
    sgn = -1;
else
    sgn = 1;
end

useNeumanProj = varCol.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end

[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotGP = 0;
if plotGP
    close all
    for patch = 1:numel(patches)
        plotNURBS(patches{patch}.nurbs,{'resolution',[100 100]});
    end
    axis equal
    axis off
    set(gca, 'Color', 'none');
    view(-100,20)
    drawnow
    hold on
    if false
        cp = zeros(size(cp_p,1),3);
        for j = 1:size(cp_p,1)
            patch = patchIdx(j);
            cp(j,:) = evaluateNURBS(patches{patch}.nurbs, cp_p(j,:));
            plot3(cp(j,1),cp(j,2),cp(j,3), '*', 'color','red')
        end
    end
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    h = findobj('type','line');
    noLines = numel(h);
    % keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_en = (p_xi+1)*(p_eta+1);
% 
p_max = max(p_xi,p_eta);
[Q2D_2,W2D_2] = tensorQuad(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
W2D_2 = flipud(W2D_2); % to reduce round-off errors in summation?
Q2D_2 = flipud(Q2D_2); % to reduce round-off errors in summation?
if quadMethodBEMsimpson
    [Q,W] = tensorQuad(p_xi+1+extraGP,p_eta+1+extraGP);
else
    [W2D,Q2D] = tensorQuad(p_xi+1+extraGP,p_eta+1+extraGP);
    noqpMax = 200;
    W = cell(noqpMax,1);
    Q = cell(noqpMax,1);
    for ii = 1:noqpMax
        [Q{ii},W{ii}] = getQuadFromFile(ii);
%         [W2D{ii},Q2D{ii}] = gaussianQuadNURBS(p_xi+1+ii+extraGP,p_eta+1+ii+extraGP);
    end
end
% [W2D,Q2D] = gaussianQuadNURBS(2*p_xi+1,2*p_eta+1);
% p_max = max(p_xi,p_eta);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(6*p_max+1,6*p_max+1);

idxRow = zeros(n_en, noElems);
Avalues = complex(zeros(n_en, noDofs, noElems)); 
Fvalues = complex(zeros(n_en, noElems, no_angles)); 
totNoQP = 0;
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
    
    J_2_x = 0.25*(Xi_e_x(2)-Xi_e_x(1))*(Eta_e_x(2)-Eta_e_x(1));
    
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
        J_1_x = norm(crossProd);

        x = R_x*pts_x;
        fact_x = J_1_x*J_2_x*wt_x;
    
        if useHBIE
            h_xi = norm(m_1);
            h_eta = norm(m_2);
            e_xi = m_1/h_xi;
            e_eta = m_2/h_eta;

            v_1 = e_xi;
            nx = crossProd/J_1_x;
            v_3 = nx;
            v_2 = cross(v_3,v_1);
            cosT = dot(e_xi,e_eta);
            sinT = dot(v_2,e_eta);
            dXIdv = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
            dR_xdv = dXIdv*[dR_xdxi; dR_xdeta];
            v_dR_xdv = [v_1.', v_2.']*dR_xdv;
        end
        if useNeumanProj
            if SHBC
                if useCBIE
                    p_inc_x = R_x*U(sctr_x,:);
                end
                if useHBIE
                    dp_inc_x = R_x*dU(sctr_x,:);
                end
            else
                dpdn_x = R_x*U(sctr_x,:);
            end
        else
            if SHBC
                if useCBIE
                    p_inc_x = p_inc(x);
                end
                if useHBIE
                    dp_inc_x = dp_inc(x,nx);
                end
            else
                if useHBIE
                    dpdn_x = dpdn(x,nx);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         keyboard
        if plotGP
            foundMarker = true;
            while foundMarker
                foundMarker = false;
                h = findobj('type','line');
                for i_h = 1:numel(h)
                    if strcmp(h(i_h).Marker,'*')
                        delete(h(i_h));
                        foundMarker = true;
                        break
                    end
                end
            end
            plot3(x(1),x(2),x(3), '*', 'color','red')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dPhi_0dny_integral = 0;
        d2Phi_0dnxdny_integral = 0;
        ugly_integral = zeros(3,1);
        FF_temp = complex(zeros(1, no_angles));
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
                noGp = size(Q2D_2,1);
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

                    rho_t = parent2ParametricSpace([0, 1],    Q2D_2(:,2));
                    theta = parent2ParametricSpace(thetaRange,Q2D_2(:,1));
                    switch area{1}
                        case 'South'
                            rho_hat = (-1 - eta_x_t)./sin(theta);
                        case 'East'
                            rho_hat = ( 1 - xi_x_t)./cos(theta);
                        case 'North'
                            rho_hat = ( 1 - eta_x_t)./sin(theta);
                        case 'West'
                            rho_hat = (-1 - xi_x_t)./cos(theta);
                    end
                    rho = rho_hat.*rho_t;

                    xi_y_t  = xi_x_t + rho.*cos(theta);
                    eta_y_t = eta_x_t + rho.*sin(theta);

                    xi_y = parent2ParametricSpace(Xi_e_y, xi_y_t);
                    eta_y = parent2ParametricSpace(Eta_e_y, eta_y_t);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                    J1 = dR_ydxi*pts_y;
                    J2 = dR_ydeta*pts_y;
                    crossProd_y = cross(J1,J2,2);
                    J_1_y = norm2(crossProd_y);
                    ny = crossProd_y./J_1_y(:,[1,1,1]);

                    J_3_y = rho;
                    J_4_y = rho_hat;
                    J_5_y = 0.25*(thetaRange(2)-thetaRange(1));
                    fact_y = J_1_y*J_2_y.*J_3_y.*J_4_y*J_5_y.*W2D_2;

                    y = R_y*pts_y;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if plotGP
                        plot3(y(:,1),y(:,2),y(:,3),'*','color','blue')
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    xmy = x(ones(noGp,1),:)-y;
                    r = norm2(xmy);
                    
                    dPhi_0dny_ = dPhi_0dny(xmy,r,ny);
                    if useHBIE
                        dPhi_0dnx_ = dPhi_0dnx(xmy,r,nx.');
                        d2Phi_0dnxdny_ = d2Phi_0dnxdny(xmy,r,nx.',ny);
                    end
                    if ~SHBC
                        if useNeumanProj
                            dpdn_y = R_y*U(sctr,:);
                        else
                            dpdn_y = dpdn(y,ny);
                        end
                        if useCBIE
                            FF_temp = FF_temp + sum(Phi_k(r).*dpdn_y.*fact_y);
                        end
                        if useHBIE
                            FF_temp = FF_temp + alpha*sum((dPhi_kdnx(xmy,r,nx.')+dPhi_0dny_).*dpdn_y.*fact_y);
                            FF_temp = FF_temp - alpha*sum(dPhi_0dny_.*(dpdn_y-dpdn_x).*fact_y);
                            FF_temp = FF_temp - alpha*dpdn_x*sum((dPhi_0dnx_.*(ny*nx.')+dPhi_0dny_+d2Phi_0dnxdny_.*(xmy*nx.')).*fact_y);
                        end
                    end

                    if useCBIE
                        CBIE = CBIE + fact_y.'*(repmat(dPhi_kdny(xmy,r,ny),1,n_en).*R_y - dPhi_0dny_*R_x);
                    end
                    if useHBIE
                        HBIE = HBIE + ((d2Phi_kdnxdny(xmy,r,nx.',ny) - d2Phi_0dnxdny_).*fact_y).'*R_y;
                        HBIE = HBIE + (fact_y.*d2Phi_0dnxdny_).'*(R_y - R_x(ones(noGp,1),:) + xmy*v_dR_xdv);
                        HBIE = HBIE + sum((dPhi_0dnx_(:,[1,1,1]).*ny + dPhi_0dny_*nx).*fact_y(:,[1,1,1]),1)*v_dR_xdv;
                    end
                end
            else
                h = diagsMax(e_y);
                if quadMethodBEMsimpson
                    x_5 = centerPts(e_y,:);
                    l = norm(x-x_5);

                    noGp = size(Q,1);
                    n_div = round(agpBEM*h/l + 1);
                    Xi_e_y_arr  = linspace(Xi_e_y(1),Xi_e_y(2),n_div+1);
                    Eta_e_y_arr = linspace(Eta_e_y(1),Eta_e_y(2),n_div+1);
                    J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1))/n_div^2;
                    xi_y = zeros(noGp,n_div^2);
                    eta_y = zeros(noGp,n_div^2);
                    counter = 1;
                    for i_eta = 1:n_div
                        Eta_e_y_sub = Eta_e_y_arr(i_eta:i_eta+1);
                        for i_xi = 1:n_div
                            Xi_e_y_sub = Xi_e_y_arr(i_xi:i_xi+1);
                            xi_y(:,counter) = parent2ParametricSpace(Xi_e_y_sub, Q(:,1));
                            eta_y(:,counter) = parent2ParametricSpace(Eta_e_y_sub, Q(:,2));
                            counter = counter + 1;
                        end
                    end
                    xi_y = reshape(xi_y,n_div^2*noGp,1);
                    eta_y = reshape(eta_y,n_div^2*noGp,1);
                    W2D_1 = repmat(W,n_div^2,1);
                else
                    xi1  = linspace(Xi_e_y(1)+Eps,Xi_e_y(2)-Eps,10);
                    if Xi_e_y(1) < xi_x && xi_x < Xi_e_y(2)
                        xi1 = [xi1, xi_x];
                    end
                    if Xi_e_y(1) < eta_x && eta_x < Xi_e_y(2)
                        xi1 = [xi1, eta_x];
                    end
                    eta1  = linspace(Eta_e_y(1)+Eps,Eta_e_y(2)-Eps,10);
                    if Eta_e_y(1) < eta_x && eta_x < Eta_e_y(2)
                        eta1 = [eta1, eta_x];
                    end
                    if Eta_e_y(1) < xi_x && xi_x < Eta_e_y(2)
                        eta1 = [eta1, xi_x];
                    end
                    [XI1,ETA1] = meshgrid(xi1,eta1);
                    XI1 = XI1(:);
                    ETA1 = ETA1(:);
                    yy = evaluateNURBS_2ndDeriv(patches{patch_y}.nurbs, [XI1,ETA1]);
                    hh = norm2(yy-x);
                    [l, I] = min(hh);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if plotGP
                        [Xs,Ys,Zs] = sphere(noSpherePts);
                        surf(pointsRadius*Xs+yy(I,1),pointsRadius*Ys+yy(I,2),pointsRadius*Zs+yy(I,3), 'FaceColor', 'green','EdgeColor','none','LineStyle','none')
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    n_qp_xi = p_xi + 1 + round(agpBEM*h/l);
                    n_qp_eta = p_eta + 1 + round(agpBEM*h/l);
                    Q_xi = repmat(Q{n_qp_xi},n_qp_eta,1);
                    Q_eta = repmat(Q{n_qp_eta}.',n_qp_xi,1);
                    Q_eta = Q_eta(:);
                    J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));
                    xi_y = parent2ParametricSpace(Xi_e_y, Q_xi);
                    eta_y = parent2ParametricSpace(Eta_e_y, Q_eta);
                    W2D_1 = W{n_qp_xi}*W{n_qp_eta}.';
                    W2D_1 = W2D_1(:);
                end
                noGp = size(xi_y,1);

                [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                J1 = dR_ydxi*pts_y;
                J2 = dR_ydeta*pts_y;
                crossProd_y = cross(J1,J2,2);
                J_1_y = norm2(crossProd_y);
                ny = crossProd_y./J_1_y(:,[1,1,1]);

                fact_y = J_1_y*J_2_y.*W2D_1;

                y = R_y*pts_y;
                xmy = x(ones(noGp,1),:)-y;
                r = norm2(xmy);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if plotGP
                    plot3(y(:,1),y(:,2),y(:,3),'*','color','blue')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~SHBC
                    if useNeumanProj
                        dpdn_y = R_y*U(sctr_y,:);
                    else
                        dpdn_y = dpdn(y,ny);
                    end
                    if useCBIE
                        FF_temp = FF_temp + sum(Phi_k(r).*dpdn_y.*fact_y);
                    end
                    if useHBIE
                        FF_temp = FF_temp + alpha*sum(dPhi_kdnx(xmy,r,nx.').*dpdn_y.*fact_y);
                    end
                end
                dPhi_0dny_ = dPhi_0dny(xmy,r,ny);
                dPhi_0dny_integral = dPhi_0dny_integral + sum(dPhi_0dny_.*fact_y); 
                if useCBIE
                    CBIE = CBIE + (dPhi_kdny(xmy,r,ny).*fact_y).'*R_y;
                end
                if useHBIE
                    d2Phi_0dnxdny_ = d2Phi_0dnxdny(xmy,r,nx.',ny);
                    d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + sum(d2Phi_0dnxdny_.*fact_y);
                    HBIE = HBIE + (d2Phi_kdnxdny(xmy,r,nx.',ny).*fact_y).'*R_y;
                    dPhi_0dnx_ = dPhi_0dnx(xmy,r,nx.');
                    ugly_integral = ugly_integral + sum((dPhi_0dnx_(:,[1,1,1]).*ny + dPhi_0dny_*nx ...
                                                                        + d2Phi_0dnxdny_(:,[1,1,1]).*xmy).*fact_y(:,[1,1,1]),1).';
                end
            end
            idxCol(:,e_y) = sctr_y;
            if useCBIE
                A_e_temp(:,:,e_y) = A_e_temp(:,:,e_y) + R_x.'*CBIE*fact_x;
            end
            if useHBIE
                A_e_temp(:,:,e_y) = A_e_temp(:,:,e_y) + alpha*R_x.'*HBIE*fact_x;
            end
            totNoQP = totNoQP + noGp;
        end
        if useCBIE
            A_e_temp(:,:,e_x) = A_e_temp(:,:,e_x) - R_x.'*R_x*(0.5*(1-sgn) + dPhi_0dny_integral)*fact_x;
        end
        if useHBIE
            temp = (dXIdv(1,:)*(v_1*ugly_integral) + dXIdv(2,:)*(v_2*ugly_integral))*[dR_xdxi; dR_xdeta];
            A_e_temp(:,:,e_x) = A_e_temp(:,:,e_x) + alpha*R_x.'*(-R_x*d2Phi_0dnxdny_integral + temp)*fact_x;
        end
        if SHBC
            if useCBIE
                F_e = F_e - R_x.'*p_inc_x*fact_x;
            end
            if useHBIE
                F_e = F_e - alpha*R_x.'*dp_inc_x*fact_x;
            end
        else
            F_e = F_e + R_x.'*FF_temp*fact_x;
            if useHBIE
                F_e = F_e + alpha*R_x.'*dpdn_x*(dPhi_0dny_integral + 0.5*(1-sgn) - v_3*ugly_integral)*fact_x;
            end
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


varCol.totNoQP = totNoQP;
