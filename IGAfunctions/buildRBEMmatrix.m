function [A, FF] = buildRBEMmatrix(varCol)

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
% internalPts = varCol.internalPts;
dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;
model = varCol.model;
extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;

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
if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end


Phi_0 = @(r)           1./(4*pi*r);
Phi_k = @(r) exp(1i*k*r)./(4*pi*r);

dPhi_0dny = @(xmy,r,ny) Phi_0(r)./r.^2*               (xmy*ny);
dPhi_kdny = @(xmy,r,ny) Phi_k(r)./r.^2.*(1 - 1i*k*r).*(xmy*ny);

dPhi_0dnx = @(xmy,r,nx) -dPhi_0dny(xmy,r,nx);
dPhi_kdnx = @(xmy,r,nx) -dPhi_kdny(xmy,r,nx);

d2Phi_0dnxdny = @(xmy,r,nx,ny) Phi_0(r)/r^2*((nx'*ny)           - 3/r^2                *(xmy*nx)*(xmy*ny));
d2Phi_kdnxdny = @(xmy,r,nx,ny) Phi_k(r)/r^2*((nx'*ny)*(1-1i*k*r)+(k^2+3/r^2*(1i*k*r-1))*(xmy*nx)*(xmy*ny));

radialPulsation = strcmp(varCol.applyLoad, 'radialPulsation');
if radialPulsation
    no_angles = 1;
    p_inc = NaN;
    dp_inc = NaN;
else
    no_angles = length(varCol.alpha_s);
    p_inc = varCol.p_inc;
    dp_inc = varCol.dp_inc;
end
dpdn = varCol.dpdn;

%% Create collocation points based on the Greville abscissae
switch varCol.model
    case {'SS', 'SS_P', 'S1', 'S1_P', 'S1_P2', 'S3', 'S5', 'MS', 'MS_P', 'EL'}
        if 0
            eps_greville_xi = 0;
            eps_greville_eta = 0;
        else
            eps_greville_xi = 1/(2*p_xi);
            eps_greville_eta = 1/(2*p_eta);
        end
    otherwise
        eps_greville_xi = 1/(2*p_xi);
        eps_greville_eta = 1/(2*p_eta);
end
% eps_greville_xi = 1/(2*p_xi);
% eps_greville_eta = 1/(2*p_eta);
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

%                 cp_p(counter2,:) = [xi_bar-eps_greville_xi/2, eta_bar-eps_greville_eta/2];
                cp_p(counter2,:) = [xi_bar, eta_bar];
                patchIdx(counter2) = patch;
                counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
        
%         [cg_xi, grev_xi] = CauchyGalerkin(p_xi, n_xi, Xi);
%         [cg_eta, grev_eta] = CauchyGalerkin(p_eta, n_eta, Eta);
% 
%         n_cp = noDofs - length(dofsToRemove);
% 
%         cp_p = zeros(n_cp,2);
%         counter = 1;
%         counter2 = 1;
%         for j = 1:n_eta
%             eta_bar = cg_eta(j);
%             for i = 1:n_xi
%                 if ~any(dofsToRemove == counter)
%                     xi_bar = cg_xi(i);      
%                     cp_p(counter2,:) = [xi_bar, eta_bar]; 
%                     counter2 = counter2 + 1;
%                 end
%                 counter = counter + 1;
%             end
%         end
    if p_xi == 1 && p_eta == 1
        eta_bar = eta_bar + 0.5*eps_greville_eta*(Eta(j+p_eta+1)-Eta(j));
        cp_p(end,:) = [xi_bar, eta_bar];
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
%     plotNURBS(patches{patch}.nurbs);
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
% ax = gca;               % get the current axis
% ax.Clipping = 'off';    % turn clipping off
% % keyboard
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psiType = 4;
    
% keyboard
%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

% [W2D,Q2D] = gaussianQuadNURBS(p_xi+1+mod(p_xi+1,2),p_eta+1+mod(p_eta+1,2));
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);
p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(3*p_max+1,3*p_max+1);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1,p_max+1);
% [W2D_2,Q2D_2] = gaussianQuadNURBS(63,63);
% [W2D,Q2D] = gaussianQuadNURBS(63,63);

A = zeros(n_cp, noDofs);
FF = zeros(n_cp, no_angles);
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
    
    Psi1_integral = 0;
    Psi2_integral = 0;
    Psi3_integral = 0;
    Psi4_integral = 0;
    dPsi1dny_integral = 0;
    dPsi2dny_integral = 0;
    
    A_row = zeros(1, noDofs);
    
    xi_x = cp_p(i,1);
    eta_x = cp_p(i,2);
    
    xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
    eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
    e_x = sum(noElemsPatch(1:patch-1)) + xi_idx + noElementsXi*(eta_idx-1);
    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts = weights(element2(e_x,:),:); % New

    [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, wgts);
    J_temp = [dR_xdxi; dR_xdeta]*pts_x;
    m_1 = J_temp(1,:);
    m_2 = J_temp(2,:);
    crossProd_x = cross(m_1,m_2);
    x = R_x*pts_x;
    if true
        e_xi = m_1/norm(m_1);
        e_eta = m_2/norm(m_2);
    end
    
    if (eta_x == 0 || eta_x == 1) && (strcmp(model,'SS') || strcmp(model,'SS_P') || strcmp(model,'S1') || strcmp(model,'S1_P') || strcmp(model,'S1_P2') || strcmp(model,'S3') ...
            || strcmp(model,'S5')  || strcmp(model,'MS') || strcmp(model,'MS_P') || strcmp(model,'EL'))
        nx = x.'/norm(x);
    else
        nx = crossProd_x'/norm(crossProd_x);
    end

    if 0
%         iiPts = 1;
%         jiPts = 2;
%         C1 = 0;
%         C2 = 0;
%         while abs(C1) > 1e3 || abs(C1) < 1e-3 || abs(C2) > 1e3 || abs(C2) < 1e-3
%             x1 = internalPts(iiPts,:);
%             x2 = internalPts(jiPts,:);
%             r1x = norm(x1-x);
%             r2x = norm(x2-x);
%             C1 = (1i*k*r2x-1)/r2x^2*dot(x2-x,nx) - (1i*k*r1x-1)/r1x^2*dot(x1-x,nx);
%             C2 = 1 - r2x^2*(1i*k*r1x-1)*dot(x1-x,nx)/(r1x^2*(1i*k*r2x-1)*dot(x2-x,nx));
%             jiPts = jiPts + 1;
%             if jiPts > size(internalPts,1)
%                 iiPts = iiPts+1;
%                 jiPts = iiPts+1;
%             end
%         end
    else
        x1 = x - 0.5*nx.';
        x2 = x - nx.';
        r1x = norm(x1-x);
        r2x = norm(x2-x);
        C2 = (1i*k*r2x-1)/r2x^2*dot(x2-x,nx) - (1i*k*r1x-1)/r1x^2*dot(x1-x,nx);
        C1 = 1 - r2x^2*(1i*k*r1x-1)*dot(x1-x,nx)/(r1x^2*(1i*k*r2x-1)*dot(x2-x,nx));
    end
    if abs(C2) < 1e-4 || abs(C1) < 1e-4
        keyboard
    end
    Phix1x = Phi_k(r1x);
    Phix2x = Phi_k(r2x);
    
    if abs(nx(1)) < 1/sqrt(2)
        d1 = sqrt(3)/2*cross([1;0;0],nx)/sqrt(1-nx(1)^2) - nx/2;
    else
        d1 = sqrt(3)/2*cross([0;1;0],nx)/sqrt(1-nx(2)^2) - nx/2;
    end
    d2 = d1+nx;
    
    if psiType == 4
        if abs(nx(1)) < 1/sqrt(2)
            d3 = sqrt(3)/2*[1-nx(1)^2; -nx(1)*nx(2); -nx(1)*nx(3)]/sqrt(1-nx(1)^2) - nx/2;
        else
            d3 = sqrt(3)/2*[-nx(1)*nx(2); 1-nx(2)^2; -nx(2)*nx(3)]/sqrt(1-nx(2)^2) - nx/2;
        end
        d4 = d3+nx;
        if abs(norm(d1)-1) > 10*eps || abs(norm(d2)-1) > 10*eps || abs(dot(d1,nx)+1/2) > 10*eps || abs(dot(d2,nx)-1/2) > 10*eps
            keyboard
        end
        if abs(norm(d3)-1) > 10*eps || abs(norm(d4)-1) > 10*eps || abs(dot(d3,nx)+1/2) > 10*eps || abs(dot(d4,nx)-1/2) > 10*eps
            keyboard
        end
        d1 = [1;0;0];
        d2 = [0;1;0];
        d3 = [0;0;1];
        d4 = [1;1;1]/sqrt(3);
        dn =   [dot(d1,nx),     dot(d2,nx),     dot(d3,nx),     dot(d4,nx)];
        dxi =  [dot(d1,e_xi),   dot(d2,e_xi),   dot(d3,e_xi),   dot(d4,e_xi)];
        deta = [dot(d1,e_eta),  dot(d2,e_eta),  dot(d3,e_eta),  dot(d4,e_eta)];
        bb = ([1,1,1,1; dn; dxi; deta]\diag([1,1/(1i*k),1/(1i*k),1/(1i*k)])).';
    else
        Psi3 = 0;
        Psi4 = 0;
        dPsi3dny = 0;
        dPsi4dny = 0;
    end
    xd = zeros(1,3);
%     xd(2) = 1/2;
    a = norm(x-xd);
    b = dot(x-xd, nx)/a;
        
    FF_temp = zeros(1, no_angles);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     quiver3([nx(1),0,nx(1),0,0],[nx(2),0,nx(2),0,0],[nx(3),0,nx(3),0,0],...
%             [d1(1),d2(1),d3(1),d4(1),nx(1)],[d1(2),d2(2),d3(2),d4(2),nx(2)],[d1(3),d2(3),d3(3),d4(3),nx(3)],'AutoScale','off')
% 	hold on
%     quiver3([0,0,0,0],[0,0,0,0],[0,0,0,0],...
%             [d1(1),d2(1),d3(1),d4(1)],[d1(2),d2(2),d3(2),d4(2)],[d1(3),d2(3),d3(3),d4(3)],'AutoScale','off')
% 	axis equal
%     if i == 3
%         figure(42)
%         xmy = @(x,y) repmat(x,size(y,1),1)-y;
%         r = @(y) norm2(xmy(x,y));
%         r1y = @(y) norm2(xmy(x1,y));
%         r2y = @(y) norm2(xmy(x2,y));
%         Phi_kTemp = @(y) Phi_k(r(y));
%         dPhi_kTemp = @(y,ny) Phi_k(r(y))./r(y).^2.*(1 - 1i*k*r(y)).*dot(xmy(x,y),ny.',2);
%         if psiType == 1
%             Psi2 = @(y) (Phi_k(r1y(y))/Phix1x - Phi_k(r2y)/Phix2x)/C2;
%             Psi1 = @(y) Phi_k(r1y(y))/Phix1x/C1 + (1-1/C1)*Phi_k(r2y(y))/Phix2x;
%             dPhi_kTemp1 = @(y,ny) Phi_k(r1y(y))./r1y(y).^2.*(1 - 1i*k*r1y(y)).*dot(xmy(x1,y),ny.',2);
%             dPhi_kTemp2 = @(y,ny) Phi_k(r2y(y))./r2y(y).^2.*(1 - 1i*k*r2y(y)).*dot(xmy(x2,y),ny.',2);
%             dPsi2dny = @(y,ny) (dPhi_kTemp1(y,ny)/Phix1x - dPhi_kTemp2(y,ny)/Phix2x)/C2;
%             dPsi1dny = @(y,ny) dPhi_kTemp1(y,ny)/Phix1x/C1 + (1-1/C1)*dPhi_kTemp2(y,ny)/Phix2x;
%             Psi3 = @(y) zeros(size(y,1));
%             Psi4 = @(y) zeros(size(y,1));
%         elseif psiType == 2
%             exp1 = @(y) exp(1i*k*dot3(-xmy(x,y),d1));
%             exp2 = @(y) exp(1i*k*dot3(-xmy(x,y),d2));
%             Psi2 = @(y) 1i*(exp1(y)-exp2(y))/k;
%             Psi1 = @(y) (exp1(y)+exp2(y))/2;
%             dPsi2dny = @(y,ny) dot3(ny.',d2).*exp2(y) - dot3(ny.',d1).*exp1(y);
%             dPsi1dny = @(y,ny) 1i*k*(dot3(ny.',d1).*exp1(y)+dot3(ny.',d2).*exp2(y))/2;
%             Psi3 = @(y) zeros(size(y,1));
%             Psi4 = @(y) zeros(size(y,1));
%         elseif psiType == 3    
%             rd = @(y) norm2(-xmy(xd,y));    
%             Psi1 = @(y) a*cos(k*(rd(y)-a))./rd(y) + sin(k*(rd(y)-a))./(k*rd(y)); % = g
%             dPsi1dny = @(y,ny) (-a*k*sin(k*(rd(y)-a))./rd(y) - a*cos(k*(rd(y)-a))./rd(y).^2 + cos(k*(rd(y)-a))./rd(y) - sin(k*(rd(y)-a))./(k*rd(y).^2)).*dot(-xmy(xd,y),ny.',2)./rd(y);
%             Psi3 = @(y) zeros(size(y,1));
%             Psi4 = @(y) zeros(size(y,1));
%         elseif psiType == 4
%             d1 = [1;0;0];
%             d2 = [0;1;0];
%             d3 = [0;0;1];
%             d4 = [1;1;1]/sqrt(3);
%             dn =   [dot(d1,nx),     dot(d2,nx),     dot(d3,nx),     dot(d4,nx)];
%             dxi =  [dot(d1,e_xi),   dot(d2,e_xi),   dot(d3,e_xi),   dot(d4,e_xi)];
%             deta = [dot(d1,e_eta),  dot(d2,e_eta),  dot(d3,e_eta),  dot(d4,e_eta)];
%             bb = ([1,1,1,1; dn; dxi; deta]\diag([1,1/(1i*k),1/(1i*k),1/(1i*k)])).';
% %             bb = ([1,1; dn(1:2)]\diag([1,1/(1i*k)])).';
%             Psi1 = @(y) bb(1,1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(1,2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(1,3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(1,4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             Psi2 = @(y) bb(2,1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(2,2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(2,3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(2,4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             Psi3 = @(y) bb(3,1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(3,2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(3,3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(3,4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             Psi4 = @(y) bb(4,1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(4,2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(4,3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(4,4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             dPsi1dny = @(y,ny) bb(1,1)*1i*k*dot3(ny.',d1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(1,2)*1i*k*dot3(ny.',d2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(1,3)*1i*k*dot3(ny.',d3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(1,4)*1i*k*dot3(ny.',d4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             dPsi2dny = @(y,ny) bb(2,1)*1i*k*dot3(ny.',d1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(2,2)*1i*k*dot3(ny.',d2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(2,3)*1i*k*dot3(ny.',d3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(2,4)*1i*k*dot3(ny.',d4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             dPsi3dny = @(y,ny) bb(3,1)*1i*k*dot3(ny.',d1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(3,2)*1i*k*dot3(ny.',d2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(3,3)*1i*k*dot3(ny.',d3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(3,4)*1i*k*dot3(ny.',d4)*exp(1i*k*dot3(-xmy(x,y),d4));
%             dPsi4dny = @(y,ny) bb(4,1)*1i*k*dot3(ny.',d1)*exp(1i*k*dot3(-xmy(x,y),d1)) + bb(4,2)*1i*k*dot3(ny.',d2)*exp(1i*k*dot3(-xmy(x,y),d2)) + bb(4,3)*1i*k*dot3(ny.',d3)*exp(1i*k*dot3(-xmy(x,y),d3)) + bb(4,4)*1i*k*dot3(ny.',d4)*exp(1i*k*dot3(-xmy(x,y),d4));
%         end
%         p_tot = @(y) varCol.analytic(y)+varCol.p_inc(y);
%         gp_tot = varCol.gAnalytic(x)+varCol.gp_inc(x);
%         integrand = @(y,ny) (p_tot(y) - p_tot(x)*Psi1(y) - dpdn(x,nx)*Psi2(y) - (gp_tot*e_xi.')*Psi3(y) - (gp_tot*e_eta.')*Psi4(y)).*dPhi_kTemp(y,ny);
%         colorFun = @(y,ny) real(integrand(y,ny));
% %         colorFun = @(y,ny) real(varCol.analytic(y).*dPhi_kTemp(y,ny));
% %         colorFun = @(y,ny) real(Phi_k(r(y)));
%         for patch = 1:numel(patches)
%             plotNURBS(patches{patch}.nurbs,{'resolution',[100 100], 'colorFun',colorFun});
%         end
%         axis equal
%         axis off
%         set(gca, 'Color', 'none');
%         view(-100,20)
%         drawnow
%         hold on
%         cp = zeros(size(cp_p,1),3);
%         for ii = 1:size(cp_p,1)
%             patch = patchIdx(ii);
%             cp(ii,:) = evaluateNURBS(patches{patch}.nurbs, cp_p(ii,:));
%             plot3(cp(ii,1),cp(ii,2),cp(ii,3), '*', 'color','red')
%         end
%         ax = gca;               % get the current axis
%         ax.Clipping = 'off';    % turn clipping off
%         
%         figure(41)
%         noPts = 10000;
%         theta = linspace(3*pi/4,pi,noPts);
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
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     if i == 1
%                         plot3(y(1),y(2),y(3), '*', 'color','blue')
%                         keyboard
%                     end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if useEnrichedBfuns
                        R_y = R_y*exp(1i*k*dot(d_vec, y));
                    end
                    xmy = x-y;
                    x1my = x1-y;
                    x2my = x2-y;
                    r = norm(xmy);

                    r1y = norm(x1my);
                    r2y = norm(x2my);
                    if useCBIE
                        Phi_kTemp = Phi_k(r);
                    end
                    if radialPulsation
                        if useCBIE
                            FF_temp = FF_temp + Phi_kTemp*dpdn(y,ny)*fact;
                        end
                    end
                    if useCBIE
                        if psiType == 1
                            Psi2 = (Phi_k(r1y)/Phix1x - Phi_k(r2y)/Phix2x)/C2; % Psi2(x) = 0
                            Psi1 = Phi_k(r1y)/Phix1x/C1 + (1-1/C1)*Phi_k(r2y)/Phix2x; % Psi1(x) = 1
                            dPsi2dny = (dPhi_kdny(x1my,r1y,ny)/Phix1x - dPhi_kdny(x2my,r2y,ny)/Phix2x)/C2; % dPsi2dny(x) = 1
                            dPsi1dny = dPhi_kdny(x1my,r1y,ny)/Phix1x/C1 + (1-1/C1)*dPhi_kdny(x2my,r2y,ny)/Phix2x; % dPsi1dny(x) = 0
                        elseif psiType == 2
                            exp1 = exp(1i*k*dot(d1,y-x));
                            exp2 = exp(1i*k*dot(d2,y-x));
                            Psi2 = 1i*(exp1-exp2)/k;
                            Psi1 = (exp1+exp2)/2;
                            dPsi2dny = dot(d2,ny)*exp2 - dot(d1,ny)*exp1;
                            dPsi1dny = 1i*k*(dot(d1,ny)*exp1+dot(d2,ny)*exp2)/2;
                        elseif psiType == 3
                            rd = norm(y-xd);
                            Psi2 = sin(k*(rd-a))/(b*k*rd); % = f
                            Psi1 = a*cos(k*(rd-a))/rd + sin(k*(rd-a))/(k*rd);
                            dPsi2dny = a/(b*k)*(k*cos(k*(rd-a))/rd - sin(k*(rd-a))/rd^2)*dot(y-xd,nx)/rd;
                            dPsi1dny = (-a*k*sin(k*(rd-a))/rd - a*cos(k*(rd-a))/rd^2 + cos(k*(rd-a))/rd - sin(k*(rd-a))/(k*rd^2))*dot(y-xd,nx)/rd;
                        else
                            Psi1 = bb(1,1)*exp(1i*k*dot(-xmy,d1)) + bb(1,2)*exp(1i*k*dot(-xmy,d2)) + bb(1,3)*exp(1i*k*dot(-xmy,d3)) + bb(1,4)*exp(1i*k*dot(-xmy,d4));
                            Psi2 = bb(2,1)*exp(1i*k*dot(-xmy,d1)) + bb(2,2)*exp(1i*k*dot(-xmy,d2)) + bb(2,3)*exp(1i*k*dot(-xmy,d3)) + bb(2,4)*exp(1i*k*dot(-xmy,d4));
                            Psi3 = bb(3,1)*exp(1i*k*dot(-xmy,d1)) + bb(3,2)*exp(1i*k*dot(-xmy,d2)) + bb(3,3)*exp(1i*k*dot(-xmy,d3)) + bb(3,4)*exp(1i*k*dot(-xmy,d4));
                            Psi4 = bb(4,1)*exp(1i*k*dot(-xmy,d1)) + bb(4,2)*exp(1i*k*dot(-xmy,d2)) + bb(4,3)*exp(1i*k*dot(-xmy,d3)) + bb(4,4)*exp(1i*k*dot(-xmy,d4));
                            dPsi1dny = bb(1,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(1,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(1,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(1,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                            dPsi2dny = bb(2,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(2,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(2,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(2,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                            dPsi3dny = bb(3,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(3,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(3,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(3,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                            dPsi4dny = bb(4,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(4,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(4,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(4,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                        end
                        dPhi_kTemp = dPhi_kdny(xmy,r,ny);

                        Psi1_integral     = Psi1_integral    + Psi1*dPhi_kTemp*fact; 
                        Psi2_integral     = Psi2_integral    + Psi2*dPhi_kTemp*fact; 
                        Psi3_integral     = Psi3_integral    + (dPsi3dny*Phi_kTemp-Psi3*dPhi_kTemp)*fact; 
                        Psi4_integral     = Psi4_integral    + (dPsi4dny*Phi_kTemp-Psi4*dPhi_kTemp)*fact; 
                        dPsi1dny_integral = dPsi1dny_integral + dPsi1dny*Phi_kTemp*fact;
                        dPsi2dny_integral = dPsi2dny_integral + dPsi2dny*Phi_kTemp*fact;

                        CBIE = CBIE + dPhi_kTemp*R_y*fact;
                    end 
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
            n_div = round(2*h/l + 1);
            
%             if e_x == e
%                 n_div = 2;
%             else
%                 n_div = 1;
%             end
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
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         if i == 1
%                             plot3(y(1),y(2),y(3), '*', 'color','blue')
% %                             keyboard
%                         end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if useEnrichedBfuns
                            R_y = R_y*exp(1i*k*dot(d_vec, y));
                        end
                        xmy = x-y;
                        x1my = x1-y;
                        x2my = x2-y;
                        r = norm(xmy);

                        r1y = norm(x1my);
                        r2y = norm(x2my);
                        if useCBIE
                            Phi_kTemp = Phi_k(r);
                        end
                        if radialPulsation
                            if useCBIE
                                FF_temp = FF_temp + Phi_kTemp*dpdn(y,ny)*fact;
                            end
                        end
                        if useCBIE
                            if psiType == 1
                                Psi2 = (Phi_k(r1y)/Phix1x - Phi_k(r2y)/Phix2x)/C2; % Psi1(x) = 0
                                Psi1 = Phi_k(r1y)/Phix1x/C1 + (1-1/C1)*Phi_k(r2y)/Phix2x; % Psi2(x) = 1
                                dPsi2dny = (dPhi_kdny(x1my,r1y,ny)/Phix1x - dPhi_kdny(x2my,r2y,ny)/Phix2x)/C2; % dPsi1dny(x) = 1
                                dPsi1dny = dPhi_kdny(x1my,r1y,ny)/Phix1x/C1 + (1-1/C1)*dPhi_kdny(x2my,r2y,ny)/Phix2x; % dPsi2dny(x) = 0
                            elseif psiType == 2
                                exp1 = exp(1i*k*dot(d1,y-x));
                                exp2 = exp(1i*k*dot(d2,y-x));
                                Psi2 = 1i*(exp1-exp2)/k;
                                Psi1 = (exp1+exp2)/2;
                                dPsi2dny = dot(d2,ny)*exp2 - dot(d1,ny)*exp1;
                                dPsi1dny = 1i*k*(dot(d1,ny)*exp1+dot(d2,ny)*exp2)/2;
                            elseif psiType == 3
                                rd = norm(y-xd);
                                Psi2 = sin(k*(rd-a))/(b*k*rd); % = f
                                Psi1 = a*cos(k*(rd-a))/rd + sin(k*(rd-a))/(k*rd);
                                dPsi2dny = a/(b*k)*(k*cos(k*(rd-a))/rd - sin(k*(rd-a))/rd^2)*dot(y-xd,ny)/rd;
                                dPsi1dny = (-a*k*sin(k*(rd-a))/rd - a*cos(k*(rd-a))/rd^2 + cos(k*(rd-a))/rd - sin(k*(rd-a))/(k*rd^2))*dot(y-xd,ny)/rd;
                            else
                                Psi1 = bb(1,1)*exp(1i*k*dot(-xmy,d1)) + bb(1,2)*exp(1i*k*dot(-xmy,d2)) + bb(1,3)*exp(1i*k*dot(-xmy,d3)) + bb(1,4)*exp(1i*k*dot(-xmy,d4));
                                Psi2 = bb(2,1)*exp(1i*k*dot(-xmy,d1)) + bb(2,2)*exp(1i*k*dot(-xmy,d2)) + bb(2,3)*exp(1i*k*dot(-xmy,d3)) + bb(2,4)*exp(1i*k*dot(-xmy,d4));
                                Psi3 = bb(3,1)*exp(1i*k*dot(-xmy,d1)) + bb(3,2)*exp(1i*k*dot(-xmy,d2)) + bb(3,3)*exp(1i*k*dot(-xmy,d3)) + bb(3,4)*exp(1i*k*dot(-xmy,d4));
                                Psi4 = bb(4,1)*exp(1i*k*dot(-xmy,d1)) + bb(4,2)*exp(1i*k*dot(-xmy,d2)) + bb(4,3)*exp(1i*k*dot(-xmy,d3)) + bb(4,4)*exp(1i*k*dot(-xmy,d4));
                                dPsi1dny = bb(1,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(1,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(1,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(1,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                                dPsi2dny = bb(2,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(2,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(2,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(2,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                                dPsi3dny = bb(3,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(3,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(3,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(3,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                                dPsi4dny = bb(4,1)*1i*k*dot(ny,d1)*exp(1i*k*dot(-xmy,d1)) + bb(4,2)*1i*k*dot(ny,d2)*exp(1i*k*dot(-xmy,d2)) + bb(4,3)*1i*k*dot(ny,d3)*exp(1i*k*dot(-xmy,d3)) + bb(4,4)*1i*k*dot(ny,d4)*exp(1i*k*dot(-xmy,d4));
                            end
                            dPhi_kTemp = dPhi_kdny(xmy,r,ny);

                            Psi1_integral    = Psi1_integral    + Psi1*dPhi_kTemp*fact; 
                            Psi2_integral    = Psi2_integral    + Psi2*dPhi_kTemp*fact; 
                            Psi3_integral    = Psi3_integral    + (dPsi3dny*Phi_kTemp-Psi3*dPhi_kTemp)*fact; 
                            Psi4_integral    = Psi4_integral    + (dPsi4dny*Phi_kTemp-Psi4*dPhi_kTemp)*fact; 
                            dPsi1dny_integral = dPsi1dny_integral + dPsi1dny*Phi_kTemp*fact;
                            dPsi2dny_integral = dPsi2dny_integral + dPsi2dny*Phi_kTemp*fact;

                            CBIE = CBIE + dPhi_kTemp*R_y*fact;
                        end
                    end
                end
            end
        end
        if useCBIE
            for j = 1:n_en
                A_row(sctr(j)) = A_row(sctr(j)) + CBIE(j);
            end
        end
    end


    if useEnrichedBfuns
        R_x = R_x*exp(1i*k*dot(d_vec, x));
    end
    if useCBIE
        for j = 1:n_en
            if psiType == 1
                A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_x(j)*(dPsi1dny_integral - Psi1_integral);
            elseif psiType == 2
                A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_x(j)*(dPsi1dny_integral - Psi1_integral - 1);
            elseif psiType == 3
                A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_x(j)*(dPsi1dny_integral - Psi1_integral - 2*pi*(1+1i/(k*a))*(1-exp(2*1i*k*a))/(4*pi));
%                 A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_x(j)*(dPsi2dny_integral - Psi2_integral);
            elseif psiType == 4
                A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_x(j)*(dPsi1dny_integral - Psi1_integral - 1) + dR_xdxi(j)*Psi3_integral + dR_xdeta(j)*Psi3_integral;
            end
        end
    end
    A(i,:) = A_row;
    if radialPulsation
        if psiType == 3
            FF(i,:) = FF(i,:) + FF_temp + dpdn(x,nx)*(Psi2_integral - dPsi2dny_integral + 2*pi*1i/(k*b)*(1-exp(2*1i*k*a))/(4*pi));
%             FF(i,:) = FF(i,:) + FF_temp + dpdn(x,nx)*(Psi1_integral - dPsi1dny_integral);
        else
            FF(i,:) = FF(i,:) + FF_temp + dpdn(x,nx)*(Psi2_integral - dPsi2dny_integral);
        end
    else
%         FF(i,:) = FF(i,:) + FF_temp - dp_inc(x,nx).'*(Psi1_integral - dPsi1dy_integral);
        if useCBIE
            FF(i,:) = FF(i,:) - p_inc(x).';
        end
    end
%     dPhi_0dny_integral+0.5
%     d2Phi_0dnxdny_integral
%     Phi_k_integralExp-d2Phi_kdnxdny_integral
%     
%     R_o = 1;
%     errorInTotArea = abs((totArea-4*pi*R_o^2)/(4*pi*R_o^2))
end


