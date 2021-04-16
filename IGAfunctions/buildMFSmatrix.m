function task = buildMFSmatrix(task)


patches = task.varCol{1}.patches;
index = task.varCol{1}.index;
noElems = task.varCol{1}.noElems;
elRangeXi = task.varCol{1}.elRange{1};
elRangeEta = task.varCol{1}.elRange{2};
pIndex = task.varCol{1}.pIndex;

k = task.misc.omega/task.varCol{1}.c_f;

SHBC = strcmp(task.misc.BC, 'SHBC');
if SHBC
    no_angles = length(task.ffp.alpha_s);
    p_inc = task.p_inc_;
    dp_inc = task.dp_inc_;
else
    no_angles = 1;
    p_inc = NaN;
    dp_inc = NaN;
end
dpdn = task.dpdn_;
% n_cp = noDofs - length(dofsToRemove);

nProgressStepSize = ceil(noElems/1000);
progressBars = task.misc.progressBars;
if progressBars
    try
        ppm = ParforProgMon('Building MFS matrix: ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
else
    ppm = NaN;
end

if false
    p_xi = task.varCol{1}.degree(1); % assume p_xi is equal in all patches
    p_eta = task.varCol{1}.degree(2); % assume p_eta is equal in all patches
    element = task.varCol{1}.element;
    element2 = task.varCol{1}.element2;
    weights = task.varCol{1}.weights;
    controlPts = task.varCol{1}.controlPts;
    knotVecs = task.varCol{1}.knotVecs;
    noCtrlPts = task.varCol{1}.noCtrlPts;
    
    xb = [min(controlPts(:,1)),max(controlPts(:,1))];
    yb = [min(controlPts(:,2)),max(controlPts(:,2))];
    zb = [min(controlPts(:,3)),max(controlPts(:,3))];


    Lx = xb(2)-xb(1);
    Ly = yb(2)-yb(1);
    Lz = zb(2)-zb(1);

    % delta = findMaxElementDiameter(task.varCol{1}.patches)/2;
    delta = 0.08;

    x = linspace(xb(1),xb(2), ceil(Lx/delta)+1);
    y = linspace(yb(1),yb(2), ceil(Ly/delta)+1);
    z = linspace(zb(1),zb(2), ceil(Lz/delta)+1);
    [x,y,z] = ndgrid(x,y,z);
    x = x(:);
    y = y(:);
    z = z(:);
    DT = delaunayTriangulation(x,y,z);
    X = DT.Points;
    % X = [-1,-1,0;
    %      1,-1,0;
    %      -1,1,0;
    %      1,1,0;
    %      -1,-1,1;
    %      1,-1,1;
    %      -1,1,1;
    %      1,1,1];
    min_d_Xon = task.mfs.parm; % minimal distance from scatterer to points in X_exterior
    if 0
        extraPts = 4; % extra knots in mesh for plotting on scatterer


        [faces,X_on,task.varCol{1}.patches] = triangulateNURBSsurface(task.varCol{1}.patches,extraPts);

        in = intriangulation(X_on,faces,X,0);
        y_s = X(in,:);

        I3 = zeros(size(y_s,1),1,'logical');
        parfor i = 1:size(y_s,1)
            I3(i) = any(norm2(repmat(y_s(i,:),size(X_on,1),1) - X_on) < min_d_Xon);
        end
        y_s(I3,:) = [];
    else
        parms = task.mfs.parms;
    %     L = parms.L;
        R_o = parms.R_o;
    %     in = or(and(sqrt(X(:,2).^2+X(:,3).^2) < R_o-min_d_Xon,and(-L <= X(:,1), X(:,1) <= 0)), ...
    %             or(and(sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2) < R_o-min_d_Xon, X(:,1) > 0), ...
    %                and(sqrt((X(:,1)+L).^2+X(:,2).^2+X(:,3).^2) < R_o-min_d_Xon, X(:,1) < -L)));
        in = sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2) < R_o-min_d_Xon;
        y_s = X(in,:);
    end
    n_sp = size(y_s,1);
    n_vec = zeros(n_sp,3);
    x_vec = zeros(n_sp,3);
    nQuadPts = ceil(sqrt(n_sp/noElems));
    noRedundantPts = nQuadPts^2*noElems - n_sp;
    dofsToRemove = round(linspace2(1,nQuadPts^2*noElems,noRedundantPts));
    Q = [copyVector(linspace2(-1,1,nQuadPts).',nQuadPts,1), copyVector(linspace2(-1,1,nQuadPts).',nQuadPts,2)];
    counter2 = 1;
    counter = 1;
    for e = 1:noElems  
        if progressBars && mod(e,nProgressStepSize) == 0
            ppm.increment();
        end
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

        for gp = 1:size(Q,1)
            pt = Q(gp,:);
            if ~any(dofsToRemove == counter)
                xi = parent2ParametricSpace(Xi_e, pt(1));
                eta = parent2ParametricSpace(Eta_e, pt(2));

                [R_y, dR_ydxi, dR_ydeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                J = pts'*[dR_ydxi' dR_ydeta'];
                crossProd = cross(J(:,1),J(:,2)); % normal vector points inwards
                J_1 = norm(crossProd);
                n_vec(counter2,:) = crossProd/J_1;
                x_vec(counter2,:) = R_y*pts;

            counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
else
    extraGP = task.misc.extraGP;
    degree = patches{1}.nurbs.degree;
    p_eta = patches{1}.nurbs.degree(2);  
    [Q, W] = gaussTensorQuad(degree+1+extraGP(1:2));
    nQuadPts = size(Q,1);
    n_vec = zeros(nQuadPts,3,noElems);
    x_vec = zeros(nQuadPts,3,noElems);
    fact = zeros(nQuadPts,noElems);
%     Q2D = [copyVector(linspace(-1,1,nQuadPts).',nQuadPts,1), copyVector(linspace(-1,1,nQuadPts).',nQuadPts,2)];
%     Q2D = [copyVector(linspace2(-1,1,nQuadPts).',nQuadPts,1), copyVector(linspace2(-1,1,nQuadPts).',nQuadPts,2)];
%     for e = 1:noElems  
    parfor e = 1:noElems  
        if progressBars && mod(e,nProgressStepSize) == 0
            ppm.increment();
        end
        patch = pIndex(e); % New
        nurbs = patches{patch}.nurbs;

        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);
        
        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        xi = parent2ParametricSpace(Xi_e, Q(:,1));
        eta = parent2ParametricSpace(Eta_e, Q(:,2));
        [y, dydxi, dydeta] = evaluateNURBS(nurbs, [xi,eta],1);

        crossProd = cross(dydxi,dydeta); % normal vector points inwards
        J_1 = norm2(crossProd);
        n_vec(:,:,e) = crossProd./J_1(:,[1,1,1]);
        x_vec(:,:,e) = y;
        
        J_1 = norm2(crossProd);
        n_vec(:,:,e) = crossProd./J_1(:,[1,1,1]);
        x_vec(:,:,e) = y;
        fact(:,e) = J_1*J_2.*W;
    end
    n_sp = nQuadPts*noElems;
    x_vec = reshape(permute(x_vec,[1,3,2]),n_sp,3);
%     [x_vec, I] = uniquetol(x_vec,1e-10,'ByRows',true);
    n_vec = reshape(permute(n_vec,[1,3,2]),n_sp,3);
    n_sp = size(x_vec,1);
%     n_vec = n_vec(I,:);
%     y_s = x_vec - n_vec*task.misc.delta;
    if task.misc.exteriorProblem
        y_s = x_vec*(1-task.mfs.delta);
    else
        y_s = x_vec*(1+task.mfs.delta);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% patches = varCol.patches;
% noPatches = varCol.noPatches;
% for patch = 1:noPatches
%     plotNURBS(patches{patch}.nurbs,{'alphaValue',0.5});
%     drawnow
% end
% axis equal
% hold on
% ax = gca;               % get the current axis
% ax.Clipping = 'off';    % turn clipping off
% axis off
% 
% plot3(y_s(:,1),y_s(:,2),y_s(:,3),'*','color','red')
% plot3(x_vec(:,1),x_vec(:,2),x_vec(:,3),'*','color','blue')
% quiver3(x_vec(:,1),x_vec(:,2),x_vec(:,3),n_vec(:,1),n_vec(:,2),n_vec(:,3))
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta = varCol.parm;
% % for i = 1:n_cp
% parfor i = 1:n_cp    
%     patch = patchIdx(i);
%     Xi = knotVecs{patch}{1}; % New
%     Eta = knotVecs{patch}{2}; % New
%     uniqueXi = unique(Xi);
%     uniqueEta = unique(Eta);
%     noElementsXi = length(uniqueXi)-1;
%     noElementsEta = length(uniqueEta)-1;
%     
%     xi_x = cp_p(i,1);
%     eta_x = cp_p(i,2);
%     
%     xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
%     eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
%     e_x = xi_idx + noElementsXi*(eta_idx-1);
%     sctr_x = element(e_x,:);
%     
%     [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi, Eta, weights);
%     pts_x = controlPts(sctr_x,:);
%     J_temp = [dR_xdxi; dR_xdeta]*pts_x;
%     m_1 = J_temp(1,:);
%     m_2 = J_temp(2,:);
%     crossProd_x = cross(m_1,m_2);
%     
%     
%     x = R_x*pts_x;
%     x_vec(i,:) = x;
%     % a minus sign is included to have the normal vector pointing inwards
%     if (eta_x == 0 || eta_x == 1) && (strcmp(model,'SS') || strcmp(model,'SS_P') || strcmp(model,'S1') || strcmp(model,'S3') ...
%             || strcmp(model,'S5')  || strcmp(model,'MS') || strcmp(model,'MS_P') || strcmp(model,'TAP') || strcmp(model,'EL'))
%         nx = -x/norm(x);
%     else
%         nx = -crossProd_x/norm(crossProd_x);
%     end
%     n_vec(i,:) = nx;
%     y_s(i,:) = x + delta*nx;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniformly distributed source points (in a volume as opposed to a surface)
% seem not to give better results
% delta = 0.5;
% n = round((1.91*varCol.noDofs)^(1/3)/(1-delta));
% x = linspace(-1,1,n);
% [X,Y,Z] = meshgrid(x,x,x);
% y_s = [reshape(X,n^3,1),reshape(Y,n^3,1),reshape(Z,n^3,1)];
% y_s(norm2(y_s) >= 1-delta,:) = [];
% 
% n_cp = size(y_s,1);
% u = rand(n_cp,2);
% theta = acos(2*u(:,2)-1);
% phi = 2*pi*u(:,1);
% x_vec = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
% n_vec = -x_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uniformly distributed points on surface actually gives worse results!?
% n_cp = varCol.noDofs;
% u = rand(n_cp,2);
% theta = acos(2*u(:,2)-1);
% phi = 2*pi*u(:,1);
% x_vec = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
% n_vec = -x_vec;
% y_s = (1-delta)*x_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indices = zeros(size(n_sp,1),1);
% y_s2 = y_s;
% for i = 1:n_sp
%     x = x_vec(i,:);
%     xmy = x(ones(n_sp,1),:) - y_s2;
%     [~,I] = min(norm2(xmy));
%     y_s2(I,:) = inf;
%     indices(i) = I;
% end
% y_s = y_s(indices,:);


if 1
    noCoresToUse = feature('numCores');
    formulation = task.misc.formulation;
    if n_sp > 1e4 && noCoresToUse < 10
        error('This is a little bit too much?')
    end
    switch formulation
        case 'SS'
            error('Implementation of formulation=SS not finished.')
%             n_sp = 1;
%             x_vec = x_vec(1:n_sp,:);
%             n_vec = n_vec(1:n_sp,:);
            
            Eps = 10*eps;
            N = floor(sqrt(n_sp)-1);
            A = zeros(n_sp,(N+1)^2);
%             F = zeros(n_sp,no_angles);
            y = task.mfs.x_0;
            R = norm2(x_vec-y(ones(n_sp,1),:));
            theta = acos((x_vec(:,3)-y(3))./R);
            phi = atan2(x_vec(:,2)-y(2),x_vec(:,1)-y(1));
            r_hat = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
            theta_hat = [cos(theta).*cos(phi), cos(theta).*sin(phi), -sin(theta)];
            phi_hat = [-sin(phi), cos(phi), zeros(n_sp,1)];
            sinTheta = sin(theta);
            zeta = k*R;
            eta = cos(theta);
            exp_phi = exp(1i*phi*(-N:N));
            
            counter = 1;
            Prev = zeros(n_sp,1);
            for n = 0:N
                h = hankel_s(n,zeta,1);
                dh = k*(n*hankel_s(n,zeta,1)./zeta - hankel_s(n+1,zeta,1));
                P = legendre(n,eta);
                
                for m = -n:n
                    P1 = P(abs(m)+1,:).';
                    if abs(m) <= n-1
                        P2 = Prev(abs(m)+1,:).';
                    else
                        P2 = zeros(n_sp,1);
                    end
                    if n == 0
                        dP = zeros(n_sp,1);
                    else
                        dP = (eta*n.*P1 - (abs(m)+n)*P2)./(eta.^2-1);
                    end
                    if m == 0
                        temp = zeros(n_sp,1);
                    else
                        temp = 1i*m./(R.*sinTheta).*exp_phi2.*h.*P1;
                    end
                    exp_phi2 = exp_phi(:,m+N+1);
                    A(:,counter) = sum(n_vec.*(k*repmat(exp_phi2.*dh.*P1,1,3).*r_hat ...
                                               - repmat(sinTheta./R.*exp_phi2.*h.*dP,1,3).*theta_hat ...
                                               - repmat(temp,1,3).*phi_hat),2);
                    if any(A(:,counter) == 0)
                        keyboard
                    end
                    counter = counter + 1;
                end
                Prev = P;
            end
            F = -dp_inc(x_vec,n_vec);
            task.mfs.N = N;
        case 'PS'
            A = zeros(n_sp);
            F = zeros(n_sp,no_angles);
%             for i = 1:n_sp
            parfor i = 1:n_sp  
                x = x_vec(i,:);
                n = n_vec(i,:);
                xmy = x(ones(n_sp,1),:) - y_s;
                A(i,:) = dPhi_kdnx(xmy,norm2(xmy),n.',k);
                if SHBC
                    F(i,:) = -dp_inc(x,n);
                else
                    F(i,:) = dpdn(x,n);
                end
            end
    end
else
    nQuadPts = 3;
    x_vec = zeros(nQuadPts^2,3,noElems);
    n_vec = zeros(nQuadPts^2,3,noElems);
    fact = zeros(nQuadPts^2,noElems);
    [Q, W] = gaussTensorQuad([nQuadPts,nQuadPts]);  
    parfor e = 1:noElems  
%     for e = 1:noElems  
        patch = pIndex(e); % New
        nurbs = patches{patch}.nurbs;

        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);
        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        xi = parent2ParametricSpace(Xi_e, Q(:,1));
        eta = parent2ParametricSpace(Eta_e, Q(:,2));
        [y, dydxi, dydeta] = evaluateNURBS(nurbs, [xi,eta], 1);

        crossProd = cross(dydxi,dydeta); % normal vector points inwards
        J_1 = norm2(crossProd);
        n_vec(:,:,e) = crossProd./J_1(:,[1,1,1]);
        x_vec(:,:,e) = y;
        fact(:,e) = J_1*J_2.*W;
    end
    n_qp = nQuadPts^2*noElems;
    x_vec = reshape(permute(x_vec,[1,3,2]),1,n_qp,3);
    n_vec = reshape(permute(n_vec,[1,3,2]),1,n_qp,3);
    fact = reshape(fact,1,n_qp);
    
    d_vec = repmat(reshape(task.d_vec.',no_angles,1,3),1,n_qp,1);
    y_s = reshape(y_s,n_sp,1,3);
    P_inc = task.misc.P_inc;
    A = zeros(n_sp);
    F = zeros(n_sp,no_angles);
    
    xmy = repmat(x_vec,n_sp,1,1) - repmat(y_s,1,n_qp,1);
    r = sqrt(sum(xmy.^2,3));
%     dPhi_kdnx_ = task.varCol{1}.dPhi_kdny(xmy,r,
    dPhi_kdnx_ = -task.varCol{1}.Phi_k(r)./r.^2.*(1 - 1i*k*r).*sum(xmy.*repmat(n_vec,n_sp,1,1),3);
    p_inc_ = P_inc*exp(1i*sum(repmat(x_vec,no_angles,1,1).*d_vec,3)*k);
    dp_inc_ = 1i*sum(d_vec.*repmat(n_vec,no_angles,1,1),3)*k.*p_inc_;
%     for i = 1:n_sp
    parfor i = 1:n_sp  
        dPhi_kdnx_i = conj(dPhi_kdnx_(i,:));
        A(i,:) = sum(dPhi_kdnx_.*repmat(dPhi_kdnx_i.*fact,n_sp,     1,1),2);
        F(i,:) = sum(-dp_inc_  .*repmat(dPhi_kdnx_i.*fact,no_angles,1,1),2);
    end
    y_s = reshape(y_s,n_sp,3);
end
task.varCol{1}.y_s = y_s;
task.varCol{1}.A_K = A;
task.varCol{1}.FF = F;
task.varCol{1}.dofsToRemove = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% plotNURBS(task.varCol{1}.nurbs,'color',getColor(1),'resolution',[20,20],'alphaValue', 0.8);
% hold on
% x_vec = reshape(x_vec,[],3);
% y_s = reshape(y_s,[],3);
% n_vec = reshape(n_vec,[],3);
% plot3(x_vec(:,1),x_vec(:,2),x_vec(:,3),'o','color','red')
% plot3(y_s(:,1),y_s(:,2),y_s(:,3),'o','color','blue')
% quiver3(x_vec(:,1),x_vec(:,2),x_vec(:,3),n_vec(:,1),n_vec(:,2),n_vec(:,3))
% axis equal
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
