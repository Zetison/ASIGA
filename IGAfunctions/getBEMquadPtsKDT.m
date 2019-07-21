function [y,ny,fact_y] = getBEMquadPtsKDT(e_x,e_y,W2D_2,Q,W,x,xi_x_t,eta_x_t,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM)

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

if e_x == e_y % use polar integration
    thetas_m = zeros(1,4);
    diag_ne2 = (1-xi_x_t)^2 + (1-eta_x_t)^2;
    diag_se2 = (1-xi_x_t)^2 + (eta_x_t+1)^2;
    diag_nw2 = (xi_x_t+1)^2 + (1-eta_x_t)^2;
    diag_sw2 = (xi_x_t+1)^2 + (eta_x_t+1)^2;
    thetas_m(1) = (diag_ne2+diag_se2-2^2)/(2*sqrt(diag_ne2*diag_se2));
    thetas_m(2) = (diag_ne2+diag_nw2-2^2)/(2*sqrt(diag_ne2*diag_nw2));
    thetas_m(3) = (diag_sw2+diag_nw2-2^2)/(2*sqrt(diag_sw2*diag_nw2));
    thetas_m(4) = (diag_sw2+diag_se2-2^2)/(2*sqrt(diag_sw2*diag_se2));
    thetas_m(thetas_m < -1) = -1;
    thetas_m(thetas_m > 1) = 1;
    thetas_m = acos(thetas_m);
    p_max = max(p_xi,p_eta);
%             no_gpSource = p_max+1+50;
    no_gpSource = 2*(p_max+1);
    n_div = ceil(sqrt(numel(W2D_2)/4)*thetas_m/(pi/2)/no_gpSource);


    n_div_r = ceil(sqrt(numel(W2D_2)/4)/no_gpSource);
    thetas = cell(4,1);
    etas_east = linspace(eps-1,1-eps,n_div(1)+1);
    thetas{1} = atan2( etas_east-eta_x_t,  1-xi_x_t);
    xis_north = linspace(1-eps,eps-1,n_div(2)+1);
    thetas{2} = atan2( 1-eta_x_t,  xis_north-xi_x_t);
    etas_west = linspace(1-eps,eps-1,n_div(3)+1);
    thetas{3} = atan2( etas_west-eta_x_t, -1-xi_x_t);
    negativeAngles = thetas{3} < 0;
    thetas{3}(negativeAngles) = thetas{3}(negativeAngles)+2*pi;
    xis_south = linspace(eps-1,1-eps,n_div(4)+1);
    thetas{4} = atan2(-1-eta_x_t,  xis_south-xi_x_t);

    J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));
    areas = {'East', 'North', 'West','South'};


    for area = 1:numel(areas)
        switch areas{area}
            case 'East'
                if abs(xi_x_t - 1) < Eps
                    n_div(area) = 0;
                    continue
                end
%                         thetaRange(area,:) = [theta_x4 theta_x1];
            case 'North'
                if abs(eta_x_t - 1) < Eps
                    n_div(area) = 0;
                    continue
                end
%                         thetaRange(area,:) = [theta_x1 theta_x2];
            case 'West'
                if abs(xi_x_t - (-1)) < Eps
                    n_div(area) = 0;
                    continue
                end
%                         thetaRange(area,:) = [theta_x2 theta_x3+2*pi];
            case 'South'
                if abs(eta_x_t - (-1)) < Eps
                    n_div(area) = 0;
                    continue
                end
%                         thetaRange(area,:) = [theta_x3 theta_x4];
        end
    end
%     n_div_r = max(n_div);
%             n_div_r = 1;
    noSubElements = sum(n_div_r*n_div);
    xi_y = zeros(no_gpSource^2*noSubElements,1);
    eta_y = zeros(no_gpSource^2*noSubElements,1);
    J_3 = zeros(no_gpSource^2*noSubElements,1);
    J_4 = zeros(no_gpSource^2*noSubElements,1);
    J_5 = zeros(no_gpSource^2*noSubElements,1);

    Q_theta = repmat(Q{no_gpSource},no_gpSource,1);
    Q_r = repmat(Q{no_gpSource}.',no_gpSource,1);
    Q_r = Q_r(:);

    rho_t_arr  = linspace(0,1,n_div_r+1);

    counter2 = 1;
    for area = 1:4
        if n_div(area) == 0
            continue
        end
        n_div_theta = n_div(area);

        rho_t = zeros(no_gpSource^2,n_div_r*n_div_theta);
        theta = zeros(no_gpSource^2,n_div_r*n_div_theta);
        counter = 1;


        for i_theta = 1:n_div_theta
            theta_sub = thetas{area}(i_theta:i_theta+1);
            for i_rho_t = 1:n_div_r
                rho_t_sub = rho_t_arr(i_rho_t:i_rho_t+1);
                rho_t(:,counter) = parent2ParametricSpace(rho_t_sub, Q_r);
                theta(:,counter) = parent2ParametricSpace(theta_sub, Q_theta);
                counter = counter + 1;
            end
        end
        noSubElementsTri = n_div_r*n_div_theta;
        rho_t = reshape(rho_t,noSubElementsTri*no_gpSource^2,1);
        theta = reshape(theta,noSubElementsTri*no_gpSource^2,1);

        switch areas{area}
            case 'East'
                rho_hat = ( 1 - xi_x_t)./cos(theta);
            case 'North'
                rho_hat = ( 1 - eta_x_t)./sin(theta);
            case 'West'
                rho_hat = (-1 - xi_x_t)./cos(theta);
            case 'South'
                rho_hat = (-1 - eta_x_t)./sin(theta);
        end
        rho = rho_hat.*rho_t;

        xi_t  = xi_x_t + rho.*cos(theta);
        eta_t = eta_x_t + rho.*sin(theta);
        indices = counter2:counter2+numel(xi_t)-1;
        xi_y(indices) = parent2ParametricSpace(Xi_e_y, xi_t);
        eta_y(indices) = parent2ParametricSpace(Eta_e_y, eta_t);
        J_3(indices) = rho;
        J_4(indices) = rho_hat;

        temp = repmat(thetas{area}(2:end)-thetas{area}(1:end-1),n_div_r*no_gpSource^2,1);
        J_5(indices) = 0.25*reshape(temp,numel(xi_t),1)/n_div_r;

        counter2 = counter2 + numel(xi_t);
    end
    noGp = size(xi_y,1);
    W2D_2_temp = W{no_gpSource}*W{no_gpSource}.';
    W2D_2_temp = W2D_2_temp(:);
    W2D_2_1 = repmat(W2D_2_temp,noSubElements,1);
    
    [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasisVec(xi_y(1:noGp), eta_y(1:noGp), p_xi, p_eta, Xi_y, Eta_y, wgts_y);

    J1 = dR_ydxi*pts_y;
    J2 = dR_ydeta*pts_y;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    ny = crossProd./J_1(:,[1,1,1]);

    fact_y = J_1*J_2_y.*J_3(1:noGp).*J_4(1:noGp).*J_5(1:noGp).*W2D_2_1(1:noGp);
else
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
