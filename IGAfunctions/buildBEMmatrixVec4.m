function [A, FF, varCol] = buildBEMmatrixVec4(varCol)

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
model = varCol.model;
extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
agpBEM = varCol.agpBEM;
exteriorProblem = varCol.exteriorProblem;

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
if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
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

%% Create collocation points
colBEM_C0 = varCol.colBEM_C0;
if p_xi == 1 && p_eta == 1
    eps_greville_xi = 1/(2*colBEM_C0*p_xi);
    eps_greville_eta = 1/(2*colBEM_C0*p_eta);
else
    eps_greville_xi = 1/(colBEM_C0*p_xi);
    eps_greville_eta = 1/(colBEM_C0*p_eta);
end
n_cp = noDofs - length(dofsToRemove);
counter2 = 1;
counter = 1;
cp_p = zeros(n_cp,2);
patchIdx = zeros(n_cp,1);
patches = varCol.patches;
[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);
for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    
    if 1
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
useNeumanProj = varCol.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = [354,317,319]
for i = [317,319,392,354]
plotGP = 1;
plotPointsAsSpheres = 0;
pointsRadius = 6e-3;
lineColor = 'blue';
lineStyle = '-';
if plotGP
    close all
    for patch = 1:numel(patches)
        plotNURBS(patches{patch}.nurbs,{'resolution',[100 100], 'elementBasedSamples',true,'samplingDistance',0.1});
    end
    axis equal
    axis off
    set(gca, 'Color', 'none');
%     view(-100,20)
    view(10,20)
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
    camlight
    % keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eNeighbour = NaN; % to avoid transparency "bug"
createElementTopology

n_en = (p_xi+1)*(p_eta+1);

p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
W2D = cell(64,1);
Q2D = cell(64,1);
for ii = 1:64-(p_xi+1+extraGP)
    [W2D{ii},Q2D{ii}] = gaussianQuadNURBS(p_xi+1+ii+extraGP,p_eta+1+ii+extraGP);
end
W2D_2 = flipud(W2D_2); % to reduce round-off errors in summation?
Q2D_2 = flipud(Q2D_2); % to reduce round-off errors in summation?
A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));
% for i = 319
% for i = 392
% for i = 317
%     if plotGP
%         keyboard
%     end
%     if ismember(i,[97,99,116])
%         plotGP = true;
%         keyboard
%     else
%         plotGP = false;
%     end
totNoQP = 0;
% for i = 1:n_cp
% parfor i = 1:n_cp
%     totArea = 0;
    patch = patchIdx(i);
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    
    A_row = complex(zeros(1, noDofs));
    xi_x = cp_p(i,1);
    eta_x = cp_p(i,2);

    xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
    eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
    e_x = sum(noElemsPatch(1:patch-1)) + xi_idx + noElementsXi*(eta_idx-1);
    
    idXi_x = index(e_x,1);
    idEta_x = index(e_x,2);

    Xi_e_x = elRangeXi(idXi_x,:);
    Eta_e_x = elRangeEta(idEta_x,:);
    if plotGP
        if i == 354
            xi_x = parent2ParametricSpace(Xi_e_x, Q2D{1}(1,1));
            eta_x = parent2ParametricSpace(Eta_e_x, Q2D{1}(1,1));
        end
    end
    
    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts = weights(element2(e_x,:),:); % New

    if useHBIE            
        singularMapping = true;
        while singularMapping
            [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, wgts);
            x = R_x*pts_x;
            J_temp = [dR_xdxi; dR_xdeta]*pts_x;
            m_1 = J_temp(1,:);
            m_2 = J_temp(2,:);
            crossProd_x = cross(m_1,m_2);
            h_xi = norm(m_1);
            if h_xi < Eps
                eps_greville_xi = 1/(2*p_xi)*(Xi_e_x(2)-Xi_e_x(1));
                if xi_x+eps_greville_xi > Xi_e_x(2)
                    xi_x = xi_x - eps_greville_xi*(Xi_e_x(2)-Xi_e_x(1));
                else
                    xi_x = xi_x + eps_greville_xi*(Xi_e_x(2)-Xi_e_x(1));
                end
                continue
            end
            h_eta = norm(m_2);
            if h_eta < Eps
                eps_greville_eta = 1/(2*p_eta)*(Eta_e_x(2)-Eta_e_x(1));
                if eta_x+eps_greville_eta > Eta_e_x(2)
                    eta_x = eta_x - eps_greville_eta*(Eta_e_x(2)-Eta_e_x(1));
                else
                    eta_x = eta_x + eps_greville_eta*(Eta_e_x(2)-Eta_e_x(1));
                end
                continue
            end
            singularMapping = false;
        end
        e_xi = m_1/h_xi;
        e_eta = m_2/h_eta;

        v_1 = m_1/norm(m_1);
        nx = crossProd_x/norm(crossProd_x);
        v_3 = nx;
        v_2 = cross(v_3,v_1);
        cosT = dot(e_xi,e_eta);
        sinT = dot(v_2,e_eta);
        dXIdv = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
        dR_xdv = dXIdv*[dR_xdxi; dR_xdeta];
        v_dR_xdv = [v_1.', v_2.']*dR_xdv;
    else
        R_x = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, wgts);
        x = R_x*pts_x;
        nx = NaN;        
    end
%     if norm(x-[0,-1,0]) < 1e-2
% %     if norm(x-[0.44774,-0.8532,0.26754]) < 1e-1
%         keyboard
%     end
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
%     Phi_k_integralExp = 0;
%     d2Phi_kdnxdny_integral = 0;
    dPhi_0dny_integral = complex(0);
    d2Phi_0dnxdny_integral = complex(0);
    ugly_integral = complex(zeros(3,1));
    FF_temp = complex(zeros(1, no_angles));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        if plotPointsAsSpheres
            [Xs,Ys,Zs] = sphere(50);
            surf(pointsRadius*Xs+x(1),pointsRadius*Ys+x(2),pointsRadius*Zs+x(3), 'FaceColor', 'blue','EdgeColor','none','LineStyle','none')
        else
            plot3(x(1),x(2),x(3), '*', 'color','red')
        end
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [adjacentElements, xi_x_tArr,eta_x_tArr] = getAdjacentElements(e_x,xi_x,eta_x,Xi_e_x,Eta_e_x,eNeighbour,Eps);
    
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
    for e_y = 1:noElems  
        patch = pIndex(e_y); % New
        Xi = knotVecs{patch}{1}; % New
        Eta = knotVecs{patch}{2}; % New

        idXi = index(e_y,1);
        idEta = index(e_y,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);

        sctr = element(e_y,:);
        pts = controlPts(sctr,:);
        wgts = weights(element2(e_y,:)); % New      
        if useCBIE
            CBIE = complex(zeros(1, n_en));
        end
        if useHBIE
            HBIE = complex(zeros(1, n_en));
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
        [collocationPointIsInElement,idx] = ismember(e_y,adjacentElements);
        if collocationPointIsInElement % use polar integration
            noGp = size(Q2D_2,1);
            xi_x_t = xi_x_tArr(idx);
            eta_x_t = eta_x_tArr(idx);
            theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
            theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
            theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
            theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            for area = {'South', 'East', 'North', 'West'}
                switch area{1}
                    case 'South'
                        if abs(eta_x_t - (-1)) < Eps
                            continue
                        end
                        thetaRange = [theta_x3 theta_x4];
                    case 'East'
                        if abs(xi_x_t - 1) < Eps
                            continue
                        end
                        thetaRange = [theta_x4 theta_x1];
                    case 'North'
                        if abs(eta_x_t - 1) < Eps
                            continue
                        end
                        thetaRange = [theta_x1 theta_x2];
                    case 'West'
                        if abs(xi_x_t - (-1)) < Eps
                            continue
                        end
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

                xi_t  = xi_x_t + rho.*cos(theta);
                eta_t = eta_x_t + rho.*sin(theta);

                xi = parent2ParametricSpace(Xi_e, xi_t);
                eta = parent2ParametricSpace(Eta_e, eta_t);

                [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                J1 = dR_ydxi*pts;
                J2 = dR_ydeta*pts;
                crossProd = cross(J1,J2,2);
                J_1 = norm2(crossProd);
                ny = crossProd./J_1(:,[1,1,1]);

                J_3 = rho;
                J_4 = rho_hat;
                J_5 = 0.25*(thetaRange(2)-thetaRange(1));
                fact_y = J_1*J_2.*J_3.*J_4*J_5.*W2D_2;

                y = R_y*pts;
                if useEnrichedBfuns
                    temp = exp(1i*k*(y*d_vec));
                    R_y = R_y.*temp(:,ones(1,noGp));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if plotGP
                    if plotPointsAsSpheres
                        [Xs,Ys,Zs] = sphere(50);
                        for ii = 1:size(y,1)
                            surf(pointsRadius*Xs+y(ii,1),pointsRadius*Ys+y(ii,2),pointsRadius*Zs+y(ii,3), 'FaceColor', 'red','EdgeColor','none','LineStyle','none')
                        end
                    else
                        plot3(y(:,1),y(:,2),y(:,3),'*','color','blue')
                    end
                    rho_t2 = linspace2(0,1,100).';
                    theta2 = thetaRange(1);
                    switch area{1}
                        case 'South'
                            rho_hat2 = (-1 - eta_x_t)./sin(theta2);
                        case 'East'
                            rho_hat2 = ( 1 - xi_x_t)./cos(theta2);
                        case 'North'
                            rho_hat2 = ( 1 - eta_x_t)./sin(theta2);
                        case 'West'
                            rho_hat2 = (-1 - xi_x_t)./cos(theta2);
                    end
                    rho2 = rho_hat2.*rho_t2;

                    xi_t2  = xi_x_t + rho2.*cos(theta2);
                    eta_t2 = eta_x_t + rho2.*sin(theta2);
                    xi = parent2ParametricSpace(Xi_e, xi_t2);
                    eta = parent2ParametricSpace(Eta_e, eta_t2);
                    if ~(all(abs(xi-Xi_e(1)) < Eps) || all(abs(xi-Xi_e(2)) < Eps) || all(abs(eta-Eta_e(1)) < Eps) || all(abs(eta-Eta_e(2)) < Eps))
                        yy = evaluateNURBS_2ndDeriv(patches{patch}.nurbs, [xi,eta]);
                        plot3(yy(:,1),yy(:,2),yy(:,3),lineStyle,'color',lineColor)
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                xmy = x(ones(noGp,1),:)-y;
                r = norm2(xmy);
                

                dPhi_0dny_ = dPhi_0dny(xmy,r,ny);
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
%                         FF_temp = FF_temp - alpha*sum(dPhi_0dny_.*(dpdn_y-dpdn_x).*fact_y);
%                         FF_temp = FF_temp - alpha*dpdn_x*sum((dPhi_0dnx_.*(ny*nx.')+dPhi_0dny_+d2Phi_0dnxdny_.*(xmy*nx.')).*fact_y);
                    end
                end
                dPhi_0dny_ = dPhi_0dny(xmy,r,ny);
                dPhi_0dny_integral = dPhi_0dny_integral + sum(dPhi_0dny_.*fact_y); 
                if useCBIE
                    CBIE = CBIE + (dPhi_kdny(xmy,r,ny).*fact_y).'*R_y;
%                     CBIE = CBIE + fact_y.'*(repmat(dPhi_kdny(xmy,r,ny),1,n_en).*R_y - dPhi_0dny_*R_x);
                end
                if useHBIE
                    d2Phi_0dnxdny_ = d2Phi_0dnxdny(xmy,r,nx.',ny);
                    d2Phi_0dnxdny_integral = d2Phi_0dnxdny_integral + sum(d2Phi_0dnxdny_.*fact_y);
                    HBIE = HBIE + (d2Phi_kdnxdny(xmy,r,nx.',ny).*fact_y).'*R_y;
                    dPhi_0dnx_ = dPhi_0dnx(xmy,r,nx.');
                    ugly_integral = ugly_integral + sum((dPhi_0dnx_(:,[1,1,1]).*ny + dPhi_0dny_*nx ...
                                                                        + d2Phi_0dnxdny_(:,[1,1,1]).*xmy).*fact_y(:,[1,1,1]),1).';
%                     HBIE = HBIE + ((d2Phi_kdnxdny(xmy,r,nx.',ny) - d2Phi_0dnxdny_).*fact_y).'*R_y;
%                     HBIE = HBIE + (fact_y.*d2Phi_0dnxdny_).'*(R_y - R_x(ones(noGp,1),:) + xmy*v_dR_xdv);
%                     HBIE = HBIE + sum((dPhi_0dnx_(:,[1,1,1]).*ny + dPhi_0dny_*nx).*fact_y(:,[1,1,1]),1)*v_dR_xdv;
                end
            end
        else
            noGp = size(Q2D,1);
            x_5 = centerPts(e_y,:);

            l = norm(x-x_5);
            h = diagsMax(e_y);
            n_div = round(agpBEM*h/l + 1);
            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
            xi = parent2ParametricSpace(Xi_e, Q2D{n_div}(:,1));
            eta = parent2ParametricSpace(Eta_e, Q2D{n_div}(:,2));
            
            W2D_1 = W2D{n_div};
            noGp = size(xi,1);

            [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

            J1 = dR_ydxi*pts;
            J2 = dR_ydeta*pts;
            crossProd = cross(J1,J2,2);
            J_1 = norm2(crossProd);
            ny = crossProd./J_1(:,[1,1,1]);
            fact_y = J_1*J_2.*W2D_1;

            y = R_y*pts;
            if useEnrichedBfuns
                temp = exp(1i*k*(y*d_vec));
                R_y = R_y.*temp(:,ones(1,noGp));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plotGP
                if plotPointsAsSpheres
                    [Xs,Ys,Zs] = sphere(50);
                    for ii = 1:size(y,1)
                        surf(pointsRadius*Xs+y(ii,1),pointsRadius*Ys+y(ii,2),pointsRadius*Zs+y(ii,3), 'FaceColor', 'red','EdgeColor','none','LineStyle','none')
                    end
                else
                    plot3(y(:,1),y(:,2),y(:,3),'*','color','blue')
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xmy = x(ones(noGp,1),:)-y;
            r = norm2(xmy);


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
        totNoQP = totNoQP + noGp;
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
        dphidv = dXIdv*[dR_xdxi; dR_xdeta];
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
    if SHBC
        if useCBIE
            FF(i,:) = FF(i,:) - p_inc_x;
        end
        if useHBIE
            FF(i,:) = FF(i,:) - alpha*dp_inc_x;
        end
    else 
        FF(i,:) = FF(i,:) + FF_temp;
        if useHBIE
            FF(i,:) = FF(i,:) + alpha*dpdn_x*(dPhi_0dny_integral + 0.5*(1-sgn) - v_3*ugly_integral);
        end
    end
    
    if plotGP
        figureFullScreen(gcf)
        export_fig(['../../graphics/BEM/S1_4_' num2str(i) '_neqp2' num2str(extraGPBEM)], '-png', '-transparent', '-r300')
    end
end


varCol.totNoQP = totNoQP;
