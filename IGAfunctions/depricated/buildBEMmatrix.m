function [A, FF] = buildBEMmatrix(varCol)
error('Depricated. Use buildCBEMmatrix instead')

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

dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;
extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
agpBEM = varCol.agpBEM;



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
if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end


Phi_0 = @(r)           1/(4*pi*r);
Phi_k = @(r) exp(1i*k*r)/(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)/r^2*             (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)/r^2*(1 - 1i*k*r)*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

pointPulsation = strcmp(varCol.applyLoad, 'pointPulsation');
if pointPulsation
    no_angles = 1;
    p_inc = NaN;
    dpdn = varCol.dpdn;
    dp_inc = NaN;
else
    no_angles = length(varCol.alpha_s);
    p_inc = varCol.p_inc;
    dp_inc = varCol.dp_inc;
    dpdn = NaN;
end
exteriorProblem = true;
if exteriorProblem
    sgn = -1;
else
    sgn = 1;
end

%% Create collocation points based on the Greville abscissae
% switch varCol.model
%     case {'SS', 'SS_P', 'S1', 'S3', 'S5', 'MS', 'MS_P', 'EL'}
%         if false
%             eps_greville_xi = 0;
%             eps_greville_eta = 0;
%         else
%             eps_greville_xi = 1/(2*p_xi);
%             eps_greville_eta = 1/(2*p_eta);
%         end
%     otherwise
%         eps_greville_xi = 1/(2*p_xi);
%         eps_greville_eta = 1/(2*p_eta);
% end
n_cp = noDofs - length(dofsToRemove);
counter2 = 1;
counter = 1;
cp_p = zeros(n_cp,2);
patchIdx = zeros(n_cp,1);
patches = varCol.patches;
for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    
    if 1
        if p_xi == 1 && p_eta == 1
            eps_greville_xi = 1/(4*p_xi);
            eps_greville_eta = 1/(4*p_eta);
        else
            eps_greville_xi = 1/(2*p_xi);
            eps_greville_eta = 1/(2*p_eta);
        end
        for j = 1:n_eta
            eta_bar = sum(Eta(j+1:j+p_eta))/p_eta;
            if Eta(j+1) == Eta(j+p_eta)
                if Eta(j+1) == Eta(j+p_eta+1)
                    eta_bar = eta_bar - eps_greville_eta*(Eta(j+p_eta+1)-Eta(j));
                else
                    eta_bar = eta_bar + eps_greville_eta*(Eta(j+p_eta+1)-Eta(j+1));
                end
            end
            for i = 1:n_xi
                if ~any(dofsToRemove == counter)
                    xi_bar = sum(Xi(i+1:i+p_xi))/p_xi;
                    if Xi(i+1) == Xi(i+p_xi)
                        if Xi(i+1) == Xi(i+p_xi+1)
                            xi_bar = xi_bar - eps_greville_xi*(Xi(i+p_xi+1)-Xi(i));
                        else
                            xi_bar = xi_bar + eps_greville_xi*(Xi(i+p_xi+1)-Xi(i+1));
                        end
                    end

                    cp_p(counter2,:) = [xi_bar, eta_bar];
                    patchIdx(counter2) = patch;
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    else
        cg_xi = splinesGL(Xi,p_xi);
        cg_eta = splinesGL(Eta,p_eta);
%         cg_xi = CauchyGalerkin(p_xi, n_xi, Xi);
%         cg_eta = CauchyGalerkin(p_eta, n_eta, Eta);
        for j = 1:n_eta
            eta_bar = cg_eta(j);
            for i = 1:n_xi
                if ~any(dofsToRemove == counter)
                    xi_bar = cg_xi(i);
                    cp_p(counter2,:) = [xi_bar, eta_bar];
                    patchIdx(counter2) = patch;
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    end
end


% cp_p2 = cp_p;
% for i = 1:n_cp
%     patch = patchIdx(i);
%     cp_p2(i,1) = cp_p2(i,1)/4 + mod(patch-1,4)*0.25;
%     if patch > 4
%         cp_p2(i,2) = cp_p2(i,2)/2 + 0.5;
%     else
%         cp_p2(i,2) = cp_p2(i,2)/2;
%     end
% end
% cp_p2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(42)
% for patch = 1:numel(patches)
%     plotNURBS(patches{patch}.nurbs,{'resolution',[40 40]});
% end
% axis equal
% axis off
% set(gca, 'Color', 'none');
% view(-100,20)
% drawnow
% hold on
% cp = zeros(size(cp_p,1),3);
% for i = 1:size(cp_p,1)
%     patch = patchIdx(i);
%     cp(i,:) = evaluateNURBS(patches{patch}.nurbs, cp_p(i,:));
%     plot3(cp(i,1),cp(i,2),cp(i,3), '*', 'color','red')
% end
% keyboard
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_en = (p_xi+1)*(p_eta+1);


[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);
p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(3*p_max+1,3*p_max+1);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(20,20);
% [W2D,Q2D] = gaussianQuadNURBS(3*20,3*20);

A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));
% for i = 1:n_cp
parfor i = 1:n_cp
%     totArea = 0;

    patch = patchIdx(i);
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
%     Phi_k_integralExp = 0;
%     d2Phi_kdnxdny_integral = 0;
    dPhi_0dny_integral = 0;
    d2Phi_0dnxdny_integral = 0;
    ugly_integral = zeros(3,1);
    A_row = zeros(1, noDofs);

    xi_x = cp_p(i,1);
    eta_x = cp_p(i,2);

    xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
    eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
    e_x = sum(noElemsPatch(1:patch-1)) + xi_idx + noElementsXi*(eta_idx-1);
    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts = weights(element2(e_x,:),:); % New

    if useHBIE
        [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, wgts);
        J_temp = [dR_xdxi; dR_xdeta]*pts_x;
        m_1 = J_temp(1,:);
        m_2 = J_temp(2,:);
        crossProd_x = cross(m_1,m_2);
        h_xi = norm(m_1);
        h_eta = norm(m_2);
        e_xi = m_1/h_xi;
        e_eta = m_2/h_eta;

        if (eta_x == 0 || eta_x == 1) && (strcmp(model,'SS') || strcmp(model,'SS_P') || strcmp(model,'S1') || strcmp(model,'S3') ...
                || strcmp(model,'S5')  || strcmp(model,'MS') || strcmp(model,'MS_P') || strcmp(model,'EL'))
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
        R_x = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, wgts);
        nx = NaN;        
    end
    x = R_x*pts_x;
    FF_temp = zeros(1, no_angles);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             nx = x.'/norm(x);        
% % %             quiver3([nx(1),0,nx(1),0,0],[nx(2),0,nx(2),0,0],[nx(3),0,nx(3),0,0],...
% % %                     [d1(1),d2(1),d3(1),d4(1),nx(1)],[d1(2),d2(2),d3(2),d4(2),nx(2)],[d1(3),d2(3),d3(3),d4(3),nx(3)],'AutoScale','off')
% % %             hold on
% % %             quiver3([0,0,0,0],[0,0,0,0],[0,0,0,0],...
% % %                     [d1(1),d2(1),d3(1),d4(1)],[d1(2),d2(2),d3(2),d4(2)],[d1(3),d2(3),d3(3),d4(3)],'AutoScale','off')
% % %             axis equal
%     Phi_k = @(r)exp(1i*k*r)./(4*pi*r);
%     if i == 3
%         figure(42)
%         xmy = @(x,y) repmat(x,size(y,1),1)-y;
%         r = @(y) norm2(xmy(x,y));
%         r1y = @(y) norm2(xmy(x1,y));
%         r2y = @(y) norm2(xmy(x2,y));
%         Phi_kTemp = @(y) Phi_k(r(y));
%         dPhi_kTemp = @(y,ny) Phi_k(r(y))./r(y).^2.*(1 - 1i*k*r(y)).*dot(xmy(x,y),ny.',2);
%         dPhi_0Temp = @(y,ny) 1/(4*pi)./r(y).^3.*dot(xmy(x,y),ny.',2);
% 
%         p_tot = @(y) varCol.analytic(y)+varCol.p_inc(y);
%         gp_tot = varCol.gAnalytic(x)+varCol.gp_inc(x);
%         integrand = @(y,ny) p_tot(y).*dPhi_kTemp(y,ny) - p_tot(x).*dPhi_0Temp(y,ny);
%         colorFun = @(y,ny) real(integrand(y,ny));
%     %         colorFun = @(y,ny) real(varCol.analytic(y).*dPhi_kTemp(y,ny));
%     %         colorFun = @(y,ny) real(Phi_k(r(y)));
%         for patch = 1:numel(patches)
%             plotNURBS(patches{patch}.nurbs,{'resolution',[200 200], 'colorFun',colorFun});
%         end
%         axis equal
%         axis off
%         set(gca, 'Color', 'none');
%         view(-100,20)
%         drawnow
%         ax = gca;               % get the current axis
%         ax.Clipping = 'off';    % turn clipping off
%         hold on
%         cp = zeros(size(cp_p,1),3);
%         for ii = 1:size(cp_p,1)
%             patch = patchIdx(ii);
%             cp(ii,:) = evaluateNURBS(patches{patch}.nurbs, cp_p(ii,:));
%             plot3(cp(ii,1),cp(ii,2),cp(ii,3), '*', 'color','red')
%         end
% 
%         figure(41)
%         noPts = 1000;
%         theta = linspace(2*pi/4,pi,noPts);
%         phi = atan2(x(2),x(1)); 
%         X = cos(phi)*sin(theta); 
%         Y = sin(phi)*sin(theta); 
%         Z = cos(theta); 
%         plot(theta,colorFun([reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)], [reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)].'))
%         hold on
%         colorFunImag = @(y,ny) imag(integrand(y,ny));
%         plot(theta,colorFunImag([reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)], [reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)].'))
%         keyboard
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:noElems  
        patch = pIndex(e); % New
        Xi = knotVecs{patch}{1}; % New
        Eta = knotVecs{patch}{2}; % New

        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);

        sctr = element(e,:);
        pts = controlPts(sctr,:);
        wgts = weights(element2(e,:)); % New      
        if useCBIE
            CBIE = zeros(1, n_en);
        end
        if useHBIE
            HBIE = zeros(1, n_en);
        end
%         if eta_x == 0
%             xi_x = Xi_e(1);
%         elseif eta_x == 1 || (xi_x == 0 && Xi_e(2) == 1)
%             xi_x = Xi_e(2);
%         elseif xi_x == 1 && Xi_e(1) == 0
%             xi_x = 0;
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if (i == 40 || i == 42) && (e == 18 || e == 19)
%             n_plot = 1000;
%             xi_t_plot = linspace2(-1,1,n_plot);
%             eta_t_plot = xi_t_plot;
%             [Xi_t_plot,Eta_t_plot] = meshgrid(xi_t_plot,eta_t_plot);
%             dGdxi_t = zeros(n_plot,n_plot);
%             dGdeta_t = zeros(n_plot,n_plot);
% 
%             J_2_xi = 0.5*(Xi_e(2)-Xi_e(1));
%             J_2_eta = 0.5*(Eta_e(2)-Eta_e(1));
%             parfor i_xi = 1:n_plot
%                 for i_eta = 1:n_plot
%                     xi_t = Xi_t_plot(i_xi,i_eta);
%                     eta_t = Eta_t_plot(i_xi,i_eta);
%                     xi  = parent2ParametricSpace(Xi_e, xi_t);
%                     eta = parent2ParametricSpace(Eta_e,eta_t);
%                     [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
%                     J = pts'*[dR_ydxi' dR_ydeta'];
%                     y = R_y*pts;
%                     xmy = x-y;
%                     r = norm(xmy);
% %                     crossProd = cross(J(:,1),J(:,2));
% %                     J_1 = norm(crossProd);
% %                     ny = crossProd/J_1;
%                     E_xi = J(:,1);
%                     E_eta = J(:,2);
%                     dXdxi_t = E_xi*J_2_xi;
%                     dXdeta_t = E_eta*J_2_eta;
%                     grad_g = -(x-y)/r + d_vec';
%                     dgdxi_t = dot(grad_g, dXdxi_t);
%                     dgdeta_t = dot(grad_g, dXdeta_t);
%                     dGdxi_t(i_xi,i_eta) = dgdxi_t;
%                     dGdeta_t(i_xi,i_eta) = dgdeta_t;
%                 end
%             end
%             figure(2)
%             quiver(Xi_t_plot(1:50:end,1:50:end),Eta_t_plot(1:50:end,1:50:end),dGdxi_t(1:50:end,1:50:end),dGdeta_t(1:50:end,1:50:end))
%             hold on
%             plot([-1,-1],[-1,1],'black')
%             plot([1,1],[-1,1],'black')
%             plot([-1,1],[-1,-1],'black')
%             plot([-1,1],[1,1],'black')
%             xlim([-1.5,1.5])
%             ylim([-1.5,1.5])
%             savefig(['temp/i_' num2str(i) 'e_' num2str(e) '.fig'])
%             hold off
%             figure(3)
%             surf(xi_t_plot,eta_t_plot, log10(sqrt(dGdxi_t.^2 + dGdeta_t.^2)),'EdgeColor','none','LineStyle','none')
%             view(0,90)
%             colorbar
%             savefig(['temp/i_' num2str(i) 'e_' num2str(e) '_img.fig'])
% %             keyboard
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if false
        if e_x == e
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

                    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                    J = pts'*[dR_ydxi' dR_ydeta'];
                    crossProd = cross(J(:,1),J(:,2));
                    J_1 = norm(crossProd);
                    ny = crossProd/J_1;

                    J_3 = rho;
                    J_4 = rho_hat;
                    J_5 = 0.25*(thetaRange(2)-thetaRange(1));
                    fact = J_1*J_2*J_3*J_4*J_5*wt;

                    y = R_y*pts;
                    if useEnrichedBfuns
                        R_y = R_y*exp(1i*k*(y*d_vec));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     hold on
%                     plot3(X(1),X(2),X(3),'*')
%                     hold off
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xmy = x-y;
                    r = norm(xmy);


                    if pointPulsation
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
            nurbs = patches{patch}.nurbs;
            x_1 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(1)+eps]).';
            x_2 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(1)+eps]).';
            x_3 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(2)-eps]).';
            x_4 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(2)-eps]).';
            x_5 = evaluateNURBS(nurbs, [mean(Xi_e),  mean(Eta_e)]).';

            l = norm(x-x_5);
            h_1 = norm(x_1-x_3);
            h_2 = norm(x_2-x_4);
            h = max(h_1,h_2);
            n_div = round(agpBEM*h/l + 1);
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
                        [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                        J = pts'*[dR_ydxi' dR_ydeta'];
                        crossProd = cross(J(:,1),J(:,2));
                        J_1 = norm(crossProd);
                        ny = crossProd/J_1;
                        fact = J_1*J_2*wt;

                        y = R_y*pts;
                        if useEnrichedBfuns
                            R_y = R_y*exp(1i*k*(y*d_vec));
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         hold on
%                         plot3(X(1),X(2),X(3),'*')
%                         hold off
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        xmy = x-y;
                        r = norm(xmy);

                        if pointPulsation
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
    if useEnrichedBfuns
        R_x = R_x*exp(1i*k*dot(d_vec, x));
    end
    if useCBIE
        for j = 1:n_en
            A_row(sctr_x(j)) = A_row(sctr_x(j)) - R_x(j)*(0.5*(1-sgn) + dPhi_0dny_integral);
%             A_row(sctr_x(j)) = A_row(sctr_x(j)) - R_x(j)*0.5;
        end
    end
    if useHBIE
        dphidv = J_x*[dR_xdxi; dR_xdeta];
        gR_x = dphidv(1,:).'*v_1 + dphidv(2,:).'*v_2;
        if useEnrichedBfuns
            gR_x = exp(1i*k*dot(d_vec, x))*gR_x + 1i*k*(d_vec*R_x).';
        end
        temp = gR_x*ugly_integral;
%         temp = (J_x(1,:)*(v_1*ugly_integral) + J_x(2,:)*(v_2*ugly_integral))*[dR_xdxi; dR_xdeta];
        for j = 1:n_en
            A_row(sctr_x(j)) = A_row(sctr_x(j)) + alpha*(-R_x(j)*d2Phi_0dnxdny_integral + temp(j));
%             A_row(sctr_x(j)) = A_row(sctr_x(j)) - R_x(j)*dPhi_0dny_integral;
        end
    end
    A(i,:) = A_row;
    if pointPulsation
        FF(i,:) = FF(i,:) + FF_temp;
        if useHBIE
            FF(i,:) = FF(i,:) + alpha*dpdn(x,nx)*(dPhi_0dny_integral + 0.5*(1-sgn) - v_3*ugly_integral);
        end
    else
        if useCBIE
            FF(i,:) = FF(i,:) - p_inc(x).';
        end
        if useHBIE
            FF(i,:) = FF(i,:) - alpha*dp_inc(x,nx).';
        end
    end
%     dPhi_0dny_integral+0.5
%     d2Phi_0dnxdny_integral
%     Phi_k_integralExp-d2Phi_kdnxdny_integral
%     
%     R_o = 1;
%     errorInTotArea = abs((totArea-4*pi*R_o^2)/(4*pi*R_o^2))
end



