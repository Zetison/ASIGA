function [relL2Error, relH1Error, relH1sError, relEnergyError] = calcErrorVec(varColCell, U_cell, options)

noDomains = length(varColCell);
nodes = cell(noDomains,1);
factors = cell(noDomains,1);
u_hs = cell(noDomains,1);
du_hs = cell(noDomains,1);
for i = 1:noDomains
    U = U_cell{i};
    varCol = varColCell{i};
    %% Extract all needed data from varCol
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

    p_xi = varCol.degree(1); % assume p_xi is equal in all patches
    p_eta = varCol.degree(2); % assume p_eta is equal in all patches
    p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches


    noCtrlPts = varCol.noCtrlPts;

    %% Preallocation and initiallizations
%     if noElems < 600
%         noQuadPts = 16;
%     else
%         noQuadPts = 10;
%     end
%     noQuadPts = 6;
%     [W3D,Q3D] = gaussianQuadNURBS(noQuadPts,noQuadPts,noQuadPts); 
%     [W3D,Q3D] = gaussianQuadNURBS(p_xi+5,p_eta+5,p_zeta+5); 
    extraGP = varCol.extraGP;
    [W3D,Q3D] = gaussianQuadNURBS(p_xi+3+extraGP,p_eta+3+extraGP,p_zeta+3+extraGP); 

    if mod(i,2) == 0
        noComponents = 3;
        noComponentsDeriv = 9;
%         noComponentsDeriv = 6;
        
        Ux = U(1:3:3*noCtrlPts);
        Uy = U(2:3:3*noCtrlPts);
        Uz = U(3:3:3*noCtrlPts);
    else
        Ux = NaN;
        Uy = NaN;
        Uz = NaN;
        noComponents = 1;
        noComponentsDeriv = 3;
    end
    u_h = zeros(size(W3D,1),noElems,noComponents);
    du_h = zeros(size(W3D,1),noElems,noComponentsDeriv);
    fact = zeros(size(W3D,1),noElems);
    points = zeros(size(W3D,1),noElems,3);
    n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);

%     for e = 1:noElems
    parfor e = 1:noElems
        patch = pIndex(e); % New
        Xi = knotVecs{patch}{1}; % New
        Eta = knotVecs{patch}{2}; % New
        Zeta = knotVecs{patch}{3}; % New
        idXi = index(e,1);
        idEta = index(e,2);
        idZeta = index(e,3);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);
        Zeta_e = elRangeZeta(idZeta,:);

        J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

        sctr = element(e,:);
        if mod(i,2) == 0
            Uxyz = [Ux(sctr) Uy(sctr) Uz(sctr)];
        else
            U_sctr = U(sctr,:);
        end
        pts = controlPts(sctr,:);
        wgts = weights(element2(e,:),:); % New

        xi   = parent2ParametricSpace(Xi_e,  Q3D(:,1));
        eta  = parent2ParametricSpace(Eta_e, Q3D(:,2));
        zeta = parent2ParametricSpace(Zeta_e, Q3D(:,3));
        [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVec(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);
        dXdxi = dRdxi*pts;
        dXdeta = dRdeta*pts;
        dXdzeta = dRdzeta*pts;
        J_1 = dot(dXdxi,cross(dXdeta,dXdzeta,2),2);

        a11 = dXdxi(:,1);
        a21 = dXdxi(:,2);
        a31 = dXdxi(:,3);
        a12 = dXdeta(:,1);
        a22 = dXdeta(:,2);
        a32 = dXdeta(:,3);
        a13 = dXdzeta(:,1);
        a23 = dXdzeta(:,2);
        a33 = dXdzeta(:,3);
        Jinv1 = [(a22.*a33-a23.*a32)./J_1, (a23.*a31-a21.*a33)./J_1, (a21.*a32-a22.*a31)./J_1];
        Jinv2 = [(a13.*a32-a12.*a33)./J_1, (a11.*a33-a13.*a31)./J_1, (a12.*a31-a11.*a32)./J_1];
        Jinv3 = [(a12.*a23-a13.*a22)./J_1, (a13.*a21-a11.*a23)./J_1, (a11.*a22-a12.*a21)./J_1];
        dRdx = repmat(Jinv1(:,1),1,n_en).*dRdxi + repmat(Jinv1(:,2),1,n_en).*dRdeta + repmat(Jinv1(:,3),1,n_en).*dRdzeta;
        dRdy = repmat(Jinv2(:,1),1,n_en).*dRdxi + repmat(Jinv2(:,2),1,n_en).*dRdeta + repmat(Jinv2(:,3),1,n_en).*dRdzeta;
        dRdz = repmat(Jinv3(:,1),1,n_en).*dRdxi + repmat(Jinv3(:,2),1,n_en).*dRdeta + repmat(Jinv3(:,3),1,n_en).*dRdzeta;
        if mod(i,2) == 0
            u_h(:,e,:) = R*Uxyz;

            dUdx = dRdx*Uxyz;
            dUdy = dRdy*Uxyz;
            dUdz = dRdz*Uxyz;

            du_h(:,e,:) = [dUdx, dUdy, dUdz];
        else
            u_h(:,e,:) = R*U_sctr;
            du_h(:,e,:) = [dRdx*U_sctr,dRdy*U_sctr,dRdz*U_sctr];
        end
        fact(:,e) = J_1 * J_2 .* W3D;
        points(:,e,:) = R*pts;
    end
    u_hs{i} = reshape(u_h, size(u_h,1)*size(u_h,2),noComponents);
    du_hs{i} = reshape(du_h, size(du_h,1)*size(du_h,2),noComponentsDeriv);
    factors{i} = reshape(fact, size(fact,1)*size(fact,2),1);
    nodes{i} = reshape(points, size(points,1)*size(points,2),3);
end

if strcmp(varColCell{1}.applyLoad, 'radialPulsation')
    data.p = varColCell{1}.analytic(nodes{1});
    dp = varColCell{1}.gAnalytic(nodes{1});
    data.dpdx = dp(:,1);
    data.dpdy = dp(:,2);
    data.dpdz = dp(:,3);
else
    data = e3Dss(nodes,options);
end
rho_f = options.rho_f;
omega = options.omega;
c_f = options.c_f;
k = omega./c_f;

H1Error = 0;
H1sError = 0;
EnergyError = 0;
L2Error = 0;
normalizationEnergy = 0;
normalizationH1 = 0;
normalizationH1s = 0;
normalizationL2 = 0;
m = 1;
for i = 1:noDomains  
    if mod(i,2) == 0
        rho_s = options.rho_s(m);
        u = [data(m).u_x,data(m).u_y,data(m).u_z];
        du = [data(m).du_xdx,data(m).du_xdy,data(m).du_xdz,...
              data(m).du_ydx,data(m).du_ydy,data(m).du_ydz,...
              data(m).du_zdx,data(m).du_zdy,data(m).du_zdz];
        sigma = [data(m).sigma_xx,data(m).sigma_yy,data(m).sigma_zz,data(m).sigma_yz,data(m).sigma_xz,data(m).sigma_xy];
        
        C = varColCell{i}.C;
        strain_vec = (C\sigma.').';
        strain_h = [du_hs{i}(:,1),du_hs{i}(:,5),du_hs{i}(:,9),du_hs{i}(:,6)+du_hs{i}(:,8),du_hs{i}(:,7)+du_hs{i}(:,3),du_hs{i}(:,2)+du_hs{i}(:,4)];
        strain_e = strain_vec-strain_h;
        
        eCe = real(sum((strain_e*C).*conj(strain_e),2)); % the usage of real() is to remove machine epsilon imaginary part
        uCu = real(sum((strain_vec*C).*conj(strain_vec),2)); % the usage of real() is to remove machine epsilon imaginary part
                
        u_e = u-u_hs{i};
        du_e = du-du_hs{i};
        
        u2 = sum(abs(u).^2,2);
        du2 = sum(abs(du).^2,2);
        u_e2 = sum(abs(u_e).^2,2);
        du_e2 = sum(abs(du_e).^2,2);
        
        H1Error = H1Error + sum((u_e2 + du_e2).*factors{i});
        H1sError = H1sError + sum(du_e2.*factors{i});
        normalizationH1 = normalizationH1 + sum((u2 + du2).*factors{i});
        normalizationH1s = normalizationH1s + sum(du2.*factors{i});
        
        EnergyError = EnergyError + sum((eCe + rho_s*omega^2*u_e2).*factors{i});
        normalizationEnergy = normalizationEnergy + sum((uCu + rho_s*omega^2*u2).*factors{i});
        
        L2Error = L2Error + sum(u_e2.*factors{i});
        normalizationL2 = normalizationL2 + sum(u2.*factors{i});
        
        m = m + 1;
    else
        p = data(m).p;
        dp = [data(m).dpdx, data(m).dpdy, data(m).dpdz];
        p_e = p-u_hs{i};
        dp_e = dp-du_hs{i};
        
        p2 = abs(p).^2;
        dp2 = sum(abs(dp).^2,2);
        p_e2 = abs(p_e).^2;
        dp_e2 = sum(abs(dp_e).^2,2);

        H1Error = H1Error + sum((p_e2 + dp_e2).*factors{i});
        H1sError = H1sError + sum(dp_e2.*factors{i});
        normalizationH1 = normalizationH1 + sum((p2 + dp2).*factors{i});
        normalizationH1s = normalizationH1s + sum(dp2.*factors{i});

        EnergyError = EnergyError + 1/(rho_f(m)*omega^2)*sum((dp_e2 + k(m)^2*p_e2).*factors{i});
        normalizationEnergy = normalizationEnergy + 1/(rho_f(m)*omega^2)*sum((dp2 + k(m)^2*p2).*factors{i});

        L2Error = L2Error + sum(p_e2.*factors{i});
        normalizationL2 = normalizationL2 + sum(p2.*factors{i});
    end
end

relEnergyError = 100*sqrt(EnergyError/normalizationEnergy);
relH1Error = 100*sqrt(H1Error/normalizationH1);
relH1sError = 100*sqrt(H1sError/normalizationH1s);
relL2Error = 100*sqrt(L2Error/normalizationL2);
