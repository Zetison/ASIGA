function [BIE, integrals, FF_temp, sctr_y, noGp] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,SHBC,useCBIE,useHBIE,dpdn,U,...
                x,nx,xi_x_tArr,eta_x_tArr,xi_x,eta_x,adjacentElements,constants,psiType,useRegul,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEMsimpson,pD)
if nargin < 44
    pD.plotGP = false;
end
alpha = 1i/k;
patch_y = pIndex(e_y); % New
Xi_y = knotVecs{patch_y}{1}; % New
Eta_y = knotVecs{patch_y}{2}; % New

idXi = index(e_y,1);
idEta = index(e_y,2);

Xi_e_y = elRangeXi(idXi,:);
Eta_e_y = elRangeEta(idEta,:);

sctr_y = element(e_y,:);
pts_y = controlPts(sctr_y,:);
wgts_y = weights(element2(e_y,:)); % New    


[collocationPointIsInElement,idx] = ismember(e_y,adjacentElements);

noqpMax = numel(Q);
if collocationPointIsInElement % use polar integration
    noGp = size(Q2D_2,1);
    xi_x_t = xi_x_tArr(idx);
    eta_x_t = eta_x_tArr(idx);
    theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
    theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
    theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
    theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

    J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));
    xi_y = zeros(4*noGp,1);
    eta_y = zeros(4*noGp,1);
    J_3 = zeros(4*noGp,1);
    J_4 = zeros(4*noGp,1);
    J_5 = zeros(4*noGp,1);
    counter = 0;
    for area = 1:4 %{'East', 'North', 'West', 'South'}
        switch area
            case 1 %'East'
                if abs(xi_x_t - 1) < Eps
                    continue
                end
                thetaRange = [theta_x4 theta_x1];
            case 2 %'North'
                if abs(eta_x_t - 1) < Eps
                    continue
                end
                thetaRange = [theta_x1 theta_x2];
            case 3 %'West'
                if abs(xi_x_t - (-1)) < Eps
                    continue
                end
                if theta_x3 < 0
                    thetaRange = [theta_x2 theta_x3+2*pi];
                else
                    thetaRange = [theta_x2 theta_x3];
                end
            case 4 %'South'
                if abs(eta_x_t - (-1)) < Eps
                    continue
                end
                thetaRange = [theta_x3 theta_x4];
        end

        rho_t = parent2ParametricSpace([0, 1],    Q2D_2(:,1));
        theta = parent2ParametricSpace(thetaRange,Q2D_2(:,2));
        switch area
            case 1 %'East'
                rho_hat = ( 1 - xi_x_t)./cos(theta);
            case 2 %'North'
                rho_hat = ( 1 - eta_x_t)./sin(theta);
            case 3 %'West'
                rho_hat = (-1 - xi_x_t)./cos(theta);
            case 4 %'South'
                rho_hat = (-1 - eta_x_t)./sin(theta);
        end
        rho = rho_hat.*rho_t;

        xi_t  = xi_x_t + rho.*cos(theta);
        eta_t = eta_x_t + rho.*sin(theta);
        indices = counter+1:counter+noGp;
        xi_y(indices) = parent2ParametricSpace(Xi_e_y, xi_t);
        eta_y(indices) = parent2ParametricSpace(Eta_e_y, eta_t);
        J_3(indices) = rho;
        J_4(indices) = rho_hat;
        J_5(indices) = 0.25*(thetaRange(2)-thetaRange(1));
        counter = counter + noGp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pD.plotGP
            rho_t2 = linspace2(Eps,1-Eps,100).';
            theta2 = thetaRange(1);
            switch area
                case 1 %'East'
                    rho_hat2 = ( 1 - xi_x_t)./cos(theta2);
                case 2 %'North'
                    rho_hat2 = ( 1 - eta_x_t)./sin(theta2);
                case 3 %'West'
                    rho_hat2 = (-1 - xi_x_t)./cos(theta2);
                case 4 %'South'
                    rho_hat2 = (-1 - eta_x_t)./sin(theta2);
            end
            rho2 = rho_hat2.*rho_t2;

            xi_t2  = xi_x_t + rho2.*cos(theta2);
            eta_t2 = eta_x_t + rho2.*sin(theta2);
            xi2 = parent2ParametricSpace(Xi_e_y, xi_t2);
            eta2 = parent2ParametricSpace(Eta_e_y, eta_t2);
            if ~(all(abs(xi2-Xi_e_y(1)) < Eps) || all(abs(xi2-Xi_e_y(2)) < Eps) || all(abs(eta2-Eta_e_y(1)) < Eps) || all(abs(eta2-Eta_e_y(2)) < Eps))
                yy = evaluateNURBS_2ndDeriv(pD.patches{pD.patch_y}.nurbs, [xi2,eta2]);
                plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    noGp = counter;
    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi_y(1:noGp), eta_y(1:noGp), p_xi, p_eta, Xi_y, Eta_y, wgts_y);

    J1 = dR_ydxi*pts_y;
    J2 = dR_ydeta*pts_y;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    ny = crossProd./J_1(:,[1,1,1]);

    fact_y = J_1*J_2_y.*J_3(1:noGp).*J_4(1:noGp).*J_5(1:noGp).*W2D_2(1:noGp);
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pD.plotGP
            for i_xi = 2:n_div
                xi2 = Xi_e_y_arr(i_xi);
                eta2 = linspace(Eta_e_y(1)+Eps,Eta_e_y(2)-Eps,100).';
                yy = evaluateNURBS_2ndDeriv(pD.patches{pD.patch_y}.nurbs, [xi2(ones(100,1)),eta2]);
                plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor)
            end
            for i_eta = 2:n_div
                eta2 = Eta_e_y_arr(i_eta);
                xi2 = linspace(Xi_e_y(1)+Eps,Xi_e_y(2)-Eps,100).';
                yy = evaluateNURBS_2ndDeriv(pD.patches{pD.patch_y}.nurbs, [xi2,eta2(ones(100,1))]);
                plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        if pD.plotGP
            plotGP(pD,yy(I,:),'green')
        end
        
        n_qp_xi = p_xi + 1 + round(agpBEM*h/l);
        n_qp_eta = p_eta + 1 + round(agpBEM*h/l);
        if n_qp_xi > noqpMax
            warning(['Requested number of Gauss points exceeds upper limit of stored Gauss points:\n' ...
                     'n_qp_xi = ' num2str(n_qp_xi) ' > ' num2str(noqpMax) ' = noqpMax'])
            n_qp_xi = noqpMax;
        end
        if n_qp_eta > noqpMax
            warning(['Requested number of Gauss points exceeds upper limit of stored Gauss points:\n' ...
                     'n_qp_eta = ' num2str(n_qp_eta) ' > ' num2str(noqpMax) ' = noqpMax'])
            warning('Requested number of Gauss points exceeds upper limit of stored Gauss points')
            n_qp_eta = noqpMax;
        end
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
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    ny = crossProd./J_1(:,[1,1,1]);
    fact_y = J_1*J_2_y.*W2D_1;
end
y = R_y*pts_y;
if pD.plotGP
    plotGP(pD,y,'red')
end

if useEnrichedBfuns
    temp = exp(1i*k*(y*d_vec));
    R_y = R_y.*temp(:,ones(1,noGp));
end
xmy = x(ones(noGp,1),:)-y;
r = norm2(xmy);

Phi_kTemp = Phi_k(r,k);
if ~SHBC
    if useNeumanProj
        dpdn_y = R_y*U(sctr_y,:);
    else
        dpdn_y = dpdn(y,ny);
    end
    if useCBIE
        FF_temp = FF_temp + sum(Phi_kTemp.*dpdn_y.*fact_y);
    end
    if useHBIE
        FF_temp = FF_temp + alpha*sum(dPhi_kdnx(xmy,r,nx.',k).*dpdn_y.*fact_y);
    end
end
dPhi_0dny_ = dPhi_0dny(xmy,r,ny);
BIE = complex(zeros(1,size(R_y,2)));
if useCBIE
    dPhi_kTemp = dPhi_kdny(xmy,r,ny,k);
    BIE = (dPhi_kTemp.*fact_y).'*R_y;
end
if useHBIE
    d2Phi_0dnxdny_ = d2Phi_0dnxdny(xmy,r,nx.',ny);
    integrals{2} = integrals{2} + sum(d2Phi_0dnxdny_.*fact_y);
    BIE = BIE + alpha*(d2Phi_kdnxdny(xmy,r,nx.',ny,k).*fact_y).'*R_y;
    dPhi_0dnx_ = dPhi_0dnx(xmy,r,nx.');
    integrals{3} = integrals{3} + sum((dPhi_0dnx_(:,[1,1,1]).*ny + dPhi_0dny_*nx ...
                                                        + d2Phi_0dnxdny_(:,[1,1,1]).*xmy).*fact_y(:,[1,1,1]),1).';
end

if useRegul
    switch psiType
        case 1
            x1 = constants{1};
            C1 = constants{2};
            C2 = constants{3};
            ymx1 = y-x1(ones(noGp,1),:);
            R1 = norm2(ymx1);
            dR1dny = sum(ymx1.*ny,2)./R1;
            Psi2 = C1^2*sin(k*(R1-C1))./(C2*k*R1); % = f
            Psi1 = C1*cos(k*(R1-C1))./R1 + sin(k*(R1-C1))./(k*R1);
            dPsi2dny = C1^2/(C2*k)*(k*cos(k*(R1-C1))./R1 - sin(k*(R1-C1))./R1.^2).*dR1dny;
            dPsi1dny = (-C1*k*sin(k*(R1-C1))./R1 - C1*cos(k*(R1-C1))./R1.^2 + cos(k*(R1-C1))./R1 - sin(k*(R1-C1))./(k*R1.^2)).*dR1dny;
        case 2
            x1 = constants{1};
            x2 = constants{2};
            C1 = constants{3};
            C2 = constants{4};
            Phix1x = constants{5};
            Phix2x = constants{6};
            x1my = x1(ones(noGp,1),:)-y;
            x2my = x2(ones(noGp,1),:)-y;

            r1y = norm2(x1my);
            r2y = norm2(x2my);
            Psi2 = (Phi_k(r1y,k)/Phix1x - Phi_k(r2y,k)/Phix2x)/C2; % Psi2(x) = 0
            Psi1 = Phi_k(r1y,k)/Phix1x/C1 + (1-1/C1)*Phi_k(r2y,k)/Phix2x; % Psi1(x) = 1
            dPsi2dny = (dPhi_kdny(x1my,r1y,ny,k)/Phix1x - dPhi_kdny(x2my,r2y,ny,k)/Phix2x)/C2; % dPsi2dny(x) = 1
            dPsi1dny = dPhi_kdny(x1my,r1y,ny,k)/Phix1x/C1 + (1-1/C1)*dPhi_kdny(x2my,r2y,ny,k)/Phix2x; % dPsi1dny(x) = 0
        case 3
            d1 = constants{1};
            d2 = constants{2};
            exp1 = exp(-1i*k*sum(d1(ones(noGp,1),:).*xmy,2));
            exp2 = exp(-1i*k*sum(d2(ones(noGp,1),:).*xmy,2));
            Psi2 = 1i*(exp1-exp2)/k;
            Psi1 = (exp1+exp2)/2;
            dPsi2dny = sum(d2(ones(noGp,1),:).*ny,2).*exp2 - sum(d1(ones(noGp,1),:).*ny,2).*exp1;
            dPsi1dny = 1i*k*(sum(d1(ones(noGp,1),:).*ny,2).*exp1+sum(d2(ones(noGp,1),:).*ny,2).*exp2)/2;
    end

    integrals{1} = integrals{1} + sum(Psi1.*dPhi_kTemp.*fact_y); 
    integrals{2} = integrals{2} + sum(Psi2.*dPhi_kTemp.*fact_y); 
    integrals{3} = integrals{3} + sum(dPsi1dny.*Phi_kTemp.*fact_y);
    integrals{4} = integrals{4} + sum(dPsi2dny.*Phi_kTemp.*fact_y);
else
    integrals{1} = integrals{1} + sum(dPhi_0dny_.*fact_y); 
end


