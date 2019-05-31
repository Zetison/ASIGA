function [BIE, integrals, FF_temp, sctr_y, noGp, pD] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,SHBC,useCBIE,useHBIE,dpdn,U,...
                x,nx,xi_x_tArr,eta_x_tArr,xi_x,eta_x,adjacentElements,constants,psiType,useRegul,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM,pD)
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
            rho_t2 = linspace(Eps,1-Eps,100).';
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
                yy = evaluateNURBS_2ndDeriv(pD.patches{patch_y}.nurbs, [xi2,eta2]);
                pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
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
    switch quadMethodBEM
        case 'Simpson'
            x_5 = centerPts(e_y,:);
            l = norm(x-x_5);

            noGp = size(Q,1);
            h = diagsMax(e_y);
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
            W2D_1 = W2D_1*J_2_y;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if pD.plotGP
                for i_xi = 2:n_div
                    xi2 = Xi_e_y_arr(i_xi);
                    eta2 = linspace(Eta_e_y(1)+Eps,Eta_e_y(2)-Eps,100).';
                    yy = evaluateNURBS_2ndDeriv(pD.patches{patch_y}.nurbs, [xi2(ones(100,1)),eta2]);
                    pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
                end
                for i_eta = 2:n_div
                    eta2 = Eta_e_y_arr(i_eta);
                    xi2 = linspace(Xi_e_y(1)+Eps,Xi_e_y(2)-Eps,100).';
                    yy = evaluateNURBS_2ndDeriv(pD.patches{patch_y}.nurbs, [xi2,eta2(ones(100,1))]);
                    pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'New'
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
                pD = plotGP(pD,yy(I,:),'green');
            end

            h = diagsMax(e_y);
            n_qp_xi = p_xi + 1 + round(agpBEM*h/l);
            n_qp_eta = p_eta + 1 + round(agpBEM*h/l);
            if n_qp_xi > noqpMax
                warning(['Requested number of Gauss points exceeds upper limit of stored Gauss points: ' ...
                         'n_qp_xi = ' num2str(n_qp_xi) ' > ' num2str(noqpMax) ' = noqpMax'])
                n_qp_xi = noqpMax;
            end
            if n_qp_eta > noqpMax
                warning(['Requested number of Gauss points exceeds upper limit of stored Gauss points: ' ...
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
            W2D_1 = W2D_1*J_2_y;
        case 'Adaptive'
            maxArraySize = 1e4;
            xi_y = zeros(maxArraySize,1);
            eta_y = zeros(maxArraySize,1);
            W2D_1 = zeros(maxArraySize,1);
            h = diagsMax(e_y);
            x_5 = centerPts(e_y,:);
            Xi_eSub = Xi_e_y;
            Eta_eSub = Eta_e_y;
            noGp = 0;
            nurbs = patches{patch_y}.nurbs;
            indexMap = reshape(5:13,3,3).';
            XI = [];
            ETA = [];
            pts = [];
            W2D_1_temp = [];
            n_qp = 0;
            corner_i = 0;
            adaptiveQuad(0);
            xi_y = xi_y(1:noGp);
            eta_y = eta_y(1:noGp);
            W2D_1 = W2D_1(1:noGp);
    end
    noGp = size(xi_y,1);

    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

    J1 = dR_ydxi*pts_y;
    J2 = dR_ydeta*pts_y;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    ny = crossProd./J_1(:,[1,1,1]);
    fact_y = J_1.*W2D_1;
end
y = R_y*pts_y;
if pD.plotGP
    pD = plotGP(pD,y(y(:,1) < -54,:),'red');
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

function adaptiveQuad(recursionLevel)

l = norm(x-x_5);
n_div = agpBEM*h/l + 1;
if n_div < 2 % || recursionLevel > 200
    n_qp_xi = round((p_xi + 1)*n_div);
    n_qp_eta = round((p_eta + 1)*n_div);
    n_qp = n_qp_xi*n_qp_eta;
    Q_xi = repmat(Q{n_qp_xi},n_qp_eta,1);
    Q_eta = repmat(Q{n_qp_eta}.',n_qp_xi,1);
    Q_eta = Q_eta(:);
    W2D_1_temp = W{n_qp_xi}*W{n_qp_eta}.';
    W2D_1_temp = W2D_1_temp(:);

    xi_y(noGp+1:noGp+n_qp) = parent2ParametricSpace(Xi_eSub, Q_xi);
    eta_y(noGp+1:noGp+n_qp) = parent2ParametricSpace(Eta_eSub, Q_eta);
    W2D_1(noGp+1:noGp+n_qp) = W2D_1_temp*0.25*(Xi_eSub(2)-Xi_eSub(1))*(Eta_eSub(2)-Eta_eSub(1));
    noGp = noGp + n_qp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pD.plotGP
        Eps = 100*eps;
        xi2 = Xi_eSub(1);
        if abs(Xi_e_y(1) - xi2) > Eps
            eta2 = linspace(Eta_eSub(1)+Eps,Eta_eSub(2)-Eps,100).';
            yy = evaluateNURBS_2ndDeriv(nurbs, [xi2(ones(100,1)),eta2]);
            pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
        end
        eta2 = Eta_eSub(1);
        if abs(Eta_e_y(1) - eta2) > Eps
            xi2 = linspace(Xi_eSub(1)+Eps,Xi_eSub(2)-Eps,100).';
            yy = evaluateNURBS_2ndDeriv(nurbs, [xi2,eta2(ones(100,1))]);
            pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Xi_eSubArr = [Xi_eSub(1), mean(Xi_eSub), Xi_eSub(2)];
    Eta_eSubArr = [Eta_eSub(1), mean(Eta_eSub), Eta_eSub(2)];
    [XI,ETA] = meshgrid(Xi_eSubArr,Eta_eSubArr);
    pts = [mean(Xi_eSubArr(1:2)), mean(Eta_eSubArr(1:2));
           mean(Xi_eSubArr(2:3)), mean(Eta_eSubArr(1:2));
           mean(Xi_eSubArr(1:2)), mean(Eta_eSubArr(2:3));
           mean(Xi_eSubArr(2:3)), mean(Eta_eSubArr(2:3));
           XI(:),ETA(:)];
    y_temp = evaluateNURBS_2ndDeriv(nurbs, pts);
    cntr = 1;
    for j = 1:2
        for i = 1:2
            corner_i = indexMap(i:i+1,j:j+1);
            h = max([norm(y_temp(corner_i(1),:)-y_temp(corner_i(4),:)), norm(y_temp(corner_i(2),:)-y_temp(corner_i(3),:))]);
            x_5 = y_temp(cntr,:);
            Xi_eSub = Xi_eSubArr(i:i+1);
            Eta_eSub = Eta_eSubArr(j:j+1);
            adaptiveQuad(recursionLevel + 1);
            cntr = cntr + 1;
        end
    end
end
end
end
% 
% 
% function [xi,eta,w,pD] = adaptiveQuad(x,agpBEM,Q,W,h,x_5,Xi_e,Eta_e,pD,nurbs,p_xi,p_eta)
% 
% refineElementsPrev = 1;
% noElementsToRefinePrev = 1;
% counter2 = 0;
% resolved = false;
% level = 1;
% Xi_e = [Xi_e(1)+eps,Xi_e(2)-eps];
% Eta_e = [Eta_e(1)+eps,Eta_e(2)-eps];
% elRangeXi_subTemp = Xi_e;
% elRangeEta_subTemp = Eta_e;
% elRangeXi_sub = [];
% elRangeEta_sub = [];
% n_divArr = [];
% while ~resolved
%     noSubElems = 2^level;
%     E = reshape(1:noSubElems^2,noSubElems,noSubElems);
%     E1 = E(1:2:noSubElems-1,2:2:noSubElems);
%     E2 = E(1:2:noSubElems-1,1:2:noSubElems-1);
%     E3 = E(2:2:noSubElems,1:2:noSubElems-1);
%     E4 = E(2:2:noSubElems,2:2:noSubElems);
%     subElemMap = [E2(:),E3(:),E1(:),E4(:)];
%     resolved = true;
%     counter = 0;
%     refineElements = zeros(4*noElementsToRefinePrev,1);
%     for i = 1:noElementsToRefinePrev
%         e_sub = refineElementsPrev(i);
%         l = norm(x-x_5(i,:));
%         n_div = agpBEM*h(i)/l + 1;
% 
%         if n_div < 2
%             counter2 = counter2 + 1;
%             n_divArr = [n_divArr; n_div];
%             elRangeXi_sub = [elRangeXi_sub; elRangeXi_subTemp(i,:)];
%             elRangeEta_sub = [elRangeEta_sub; elRangeEta_subTemp(i,:)];
%         else
%             refineElements(counter+1:counter+4) = subElemMap(e_sub,:);
%             counter = counter+4;
%             resolved = false;
%         end
%     end
%     noElementsToRefinePrev = counter;
%     refineElementsPrev = refineElements;
%     
%     if ~resolved % compute data for next level
%         if level > 15
%             error(['x = ' num2str(x)])
%         end
%         noNodes = 2^level+1;
%         I = reshape(1:noNodes^2,noNodes,noNodes);
%         I1 = I(1:noNodes-1,2:noNodes);
%         I2 = I(1:noNodes-1,1:noNodes-1);
%         I3 = I(2:noNodes,1:noNodes-1);
%         I4 = I(2:noNodes,2:noNodes);
%         elems = refineElements(1:counter);
%         xi_t = linspace(-1,1,noNodes).';
%         eta_t = linspace(-1,1,noNodes).';
%         elRangeXi_subTemp = parent2ParametricSpace(Xi_e, [xi_t(mod(I2(elems)-1,noNodes)+1), xi_t(mod(I3(elems)-1,noNodes)+1)]);
%         elRangeEta_subTemp = parent2ParametricSpace(Eta_e, [eta_t(floor((I2(elems)-1)/noNodes)+1), eta_t(floor((I1(elems)-1)/noNodes)+1)]);
%         indices = unique([I1(elems); I2(elems); I3(elems); I4(elems)]);
%         xi = [parent2ParametricSpace(Xi_e, xi_t(mod(indices-1,noNodes)+1)); mean(elRangeXi_subTemp,2)];
%         eta = [parent2ParametricSpace(Eta_e, eta_t(floor((indices-1)/noNodes)+1)); mean(elRangeEta_subTemp,2)];
%         y = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
%         yy = zeros(noNodes^2,3);
%         yy(indices,:) = y(1:end-counter,:);
%         h = max([norm2(yy(I1(elems),:)-yy(I3(elems),:)), norm2(yy(I2(elems),:)-yy(I4(elems),:))],[],2);
%         x_5 = y(end-counter+1:end,:);
%     end
%     level = level + 1;
% end
% noSubElemsTot = counter2;
% n_divArr = n_divArr(1:noSubElemsTot);
% xi = cell(noSubElemsTot,1);
% eta = cell(noSubElemsTot,1);
% w = cell(noSubElemsTot,1);
% counter = 1;
% for i = 1:noSubElemsTot
%     n_div = n_divArr(i);
% 
%     Xi_e_y = elRangeXi_sub(i,:);
%     Eta_e_y = elRangeEta_sub(i,:);
%     n_qp_xi = round((p_xi + 1)*n_div);
%     n_qp_eta = round((p_eta + 1)*n_div);
%     
%     Q_xi = repmat(Q{n_qp_xi},n_qp_eta,1);
%     Q_eta = repmat(Q{n_qp_eta}.',n_qp_xi,1);
%     Q_eta = Q_eta(:);
%     W2D_1 = W{n_qp_xi}*W{n_qp_eta}.';
%     W2D_1 = W2D_1(:);
%             
%     xi{i} = parent2ParametricSpace(Xi_e_y, Q_xi);
%     eta{i} = parent2ParametricSpace(Eta_e_y, Q_eta);
%     w{i} = W2D_1*0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));
%     counter = counter+1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if pD.plotGP
%         Eps = 100*eps;
%         xi2 = Xi_e_y(1);
%         if abs(Xi_e(1) - xi2) > Eps
%             eta2 = linspace(Eta_e_y(1)+Eps,Eta_e_y(2)-Eps,100).';
%             yy = evaluateNURBS_2ndDeriv(nurbs, [xi2(ones(100,1)),eta2]);
%             pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
%         end
%         eta2 = Eta_e_y(1);
%         if abs(Eta_e(1) - eta2) > Eps
%             xi2 = linspace(Xi_e_y(1)+Eps,Xi_e_y(2)-Eps,100).';
%             yy = evaluateNURBS_2ndDeriv(nurbs, [xi2,eta2(ones(100,1))]);
%             pD.h(end+1) = plot3(yy(:,1),yy(:,2),yy(:,3),pD.lineStyle,'color',pD.lineColor);
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% xi = cell2mat(xi);
% eta = cell2mat(eta);
% w = cell2mat(w);

