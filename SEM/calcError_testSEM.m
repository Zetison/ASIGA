function [relL2Error, relH1Error, relH1sError, relEnergyError] = calcError_testSEM(varColCell, U_cell, options)

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
    noDofs = varCol.noDofs;

    p_xi = varCol.degree(1); % assume p_xi is equal in all patches
    p_eta = varCol.degree(2); % assume p_eta is equal in all patches
    p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches


    noCtrlPts = varCol.noCtrlPts;

    %% Preallocation and initiallizations
    if noElems < 600
        noQuadPts = 16;
    else
        noQuadPts = 10;
    end
    noQuadPts = 6;
    [W3D,Q3D] = gaussianQuadNURBS(noQuadPts,noQuadPts,noQuadPts); 

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

    nurbs = varCol.patches{1}.nurbs;
    Nxi = nurbs.number(1);
    Neta = nurbs.number(2);
    Nzeta = nurbs.number(3);
    GLL = nurbs.GLL;
    tXi = GLL{1};
    tEta = GLL{2};
    tZeta = GLL{3};
%     rho = nurbs.rho;
%     rho_xi = rho{1};
%     rho_eta = rho{2};
%     rho_zeta = rho{3};
    nxi = Nxi+3;
    neta = Neta+3;
    nzeta = Nzeta+3;
    [Wxi,Qxi] = gaussianQuadNURBS(nxi); 
    [Weta,Qeta] = gaussianQuadNURBS(neta); 
    [Wzeta,Qzeta] = gaussianQuadNURBS(nzeta); 
    
    Q3D = [copyVector2(Qxi.',neta,nzeta,1), copyVector2(Qeta.',nxi,nzeta,2), copyVector2(Qzeta.',nxi,neta,3)];
    W3D = copyVector2(Wxi.',neta,nzeta,1).*copyVector2(Weta.',nxi,nzeta,2).*copyVector2(Wzeta.',nxi,neta,3);
    
    u_h = zeros(size(W3D,1),noElems,noComponents);
    du_h = zeros(size(W3D,1),noElems,noComponentsDeriv);
    fact = zeros(size(W3D,1),noElems);
    points = zeros(size(W3D,1),noElems,3);
    parfor e = 1:noElems
%     for e = 1:noElems
        J_2 = 1;

        u_h_temp = zeros(size(W3D,1),noComponents);
        du_h_temp = zeros(size(W3D,1),noComponentsDeriv);
        fact_temp = zeros(size(W3D,1),1);
        points_temp = zeros(size(W3D,1),3);
        sctr = element(e,:);
        if mod(i,2) == 0
            Uxyz = [Ux(sctr) Uy(sctr) Uz(sctr)];
        else
            U_sctr = U(sctr,:);
        end
        pts = controlPts(sctr,:);

        for gp = 1:size(W3D,1)
            pt = Q3D(gp,:);
            wt = W3D(gp);

            xi   = pt(1);
            eta  = pt(2);
            zeta = pt(3);

            Rxi = zeros(1,Nxi);
            dRdxit = zeros(1,Nxi);
            for ii = 1:Nxi
                dRdxit(ii) = lagrangePolynomialsDeriv(xi,ii,Nxi,tXi);
                Rxi(ii) = lagrangePolynomials(xi,ii,Nxi,tXi);
            end
            Reta = zeros(1,Neta);
            dRdetat = zeros(1,Neta);
            for jj = 1:Neta
                dRdetat(jj) = lagrangePolynomialsDeriv(eta,jj,Neta,tEta);
                Reta(jj) = lagrangePolynomials(eta,jj,Neta,tEta);
            end
            Rzeta = zeros(1,Nzeta);
            dRdzetat = zeros(1,Nzeta);
            for ll = 1:Nzeta
                dRdzetat(ll) = lagrangePolynomialsDeriv(zeta,ll,Nzeta,tZeta);
                Rzeta(ll) = lagrangePolynomials(zeta,ll,Nzeta,tZeta);
            end
            R = (copyVector2(Rxi,Neta,Nzeta,1).*copyVector2(Reta,Nxi,Nzeta,2).*copyVector2(Rzeta,Nxi,Neta,3)).';
            dRdxi = (copyVector2(dRdxit,Neta,Nzeta,1).*copyVector2(Reta,Nxi,Nzeta,2).*copyVector2(Rzeta,Nxi,Neta,3)).';
            dRdeta = (copyVector2(Rxi,Neta,Nzeta,1).*copyVector2(dRdetat,Nxi,Nzeta,2).*copyVector2(Rzeta,Nxi,Neta,3)).';
            dRdzeta = (copyVector2(Rxi,Neta,Nzeta,1).*copyVector2(Reta,Nxi,Nzeta,2).*copyVector2(dRdzetat,Nxi,Neta,3)).';

            J = pts'*[dRdxi' dRdeta' dRdzeta'];
            J_1 = det(J);
            X = R*pts;
            dRdX = J'\[dRdxi; dRdeta; dRdzeta];
            if mod(i,2) == 0
                u_h_temp(gp,:) = R*Uxyz;
                
                dUdx = dRdX(1,:)*Uxyz;
                dUdy = dRdX(2,:)*Uxyz;
                dUdz = dRdX(3,:)*Uxyz;
                
                du_h_temp(gp,:) = [dUdx, dUdy, dUdz];
%                 du_h_temp(gp,:) = calculateStrainVector(dUdx, dUdy, dUdz);
            else
                u_h_temp(gp,:) = R*U_sctr;
                du_h_temp(gp,:) = dRdX*U_sctr;
            end
            fact_temp(gp) = J_1 * J_2 * wt;
            points_temp(gp,:) = X;
        end
        u_h(:,e,:) = u_h_temp;
        du_h(:,e,:) = du_h_temp;
        fact(:,e) = fact_temp;
        points(:,e,:) = points_temp;
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
