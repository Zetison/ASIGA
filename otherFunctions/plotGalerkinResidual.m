function plotGalerkinResidual(varCol,U)

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

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
analytic = varCol.analytic;
k = varCol.k;

Phi_0 = @(r)           1./(4*pi*r);
Phi_k = @(r) exp(1i*k*r)./(4*pi*r);


p_inc = varCol.p_inc;

npts = 200;

patches = varCol.patches;
noPatches = varCol.noPatches;
hold on
x_grev = [];
x_cg = [];
x_gl = [];
for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};

    xi1D = splinesGL(Xi,p_xi);
    eta1D = splinesGL(Eta,p_eta);    
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);    
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_gl = [x_gl; x];
    
    xi1D = aveknt(Xi,p_xi+1).';
    eta1D = aveknt(Eta,p_eta+1).';    
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);      
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_grev = [x_grev; x];

    xi1D = CauchyGalerkin(p_xi, numel(xi1D), Xi);
    eta1D = CauchyGalerkin(p_eta, numel(eta1D), Eta);
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);      
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_cg = [x_cg; x];
end

plot3(x_grev(:,1),x_grev(:,2),x_grev(:,3),'*','color','blue','DisplayName','Greville abscissae')
plot3(x_gl(:,1),x_gl(:,2),x_gl(:,3),'*','color','red','DisplayName','GL nodes')
plot3(x_cg(:,1),x_cg(:,2),x_cg(:,3),'*','color','green','DisplayName','Cauchy Galerkin nodes')

legend show
set(0,'DefaultLegendAutoUpdate','off')

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);
p_max = max(p_xi,p_eta);
[W2D_2,Q2D_2] = gaussianQuadNURBS(p_max+1+extraGPBEM,p_max+1+extraGPBEM);
for e_x = 1:noElems  
    patch = pIndex(e_x); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e_x,1);
    idEta = index(e_x,2);

    Xi_e_x = elRangeXi(idXi,:);
    Eta_e_x = elRangeEta(idEta,:);

    sctr = element(e_x,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e_x,:)); % New      
    U_sctr = U(sctr,:); % New      

    xi_x = copyVector(linspace(Xi_e_x(1)+10*eps,Xi_e_x(2)-10*eps,npts),npts,1);
    eta_x = copyVector(linspace(Eta_e_x(1)+10*eps,Eta_e_x(2)-10*eps,npts),npts,2);
    [x,u_x] = numericalSolEvalVectorized(xi_x,eta_x,p_xi,p_eta,Xi,Eta,wgts,pts,U_sctr);
    u_x = u_x + p_inc(x);
    f = p_inc(x) - u_x;
   
    
%     for e_y = 1:noElems  
    parfor e_y = 1:noElems  
        patch = pIndex(e_y); % New
        Xi = knotVecs{patch}{1}; % New
        Eta = knotVecs{patch}{2}; % New

        idXi = index(e_y,1);
        idEta = index(e_y,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);

        sctr = element(e_y,:);
        pts = controlPts(sctr,:);
        wgts = weights(element2(e_y,:)); 
        U_sctr = U(sctr,:); % New      
        
        if e_x == e_y
            xi_x_t = parametric2parentSpace(Xi_e_x, xi_x);
            eta_x_t = parametric2parentSpace(Eta_e_x, eta_x);
            theta_x1 = atan2( 1-eta_x_t,  1-xi_x_t);
            theta_x2 = atan2( 1-eta_x_t, -1-xi_x_t);
            theta_x3 = atan2(-1-eta_x_t, -1-xi_x_t);
            theta_x4 = atan2(-1-eta_x_t,  1-xi_x_t);

            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            for area = {'South', 'East', 'North', 'West'}
                switch area{1}
                    case 'South'
                        indices = ~(abs(eta_x - Eta_e(1)) < 100*eps);
                        thetaRange = [theta_x3 theta_x4];
                    case 'East'
                        indices = ~(abs(xi_x - Xi_e(2)) < 100*eps);
                        thetaRange = [theta_x4 theta_x1];
                    case 'North'
                        indices = ~(abs(eta_x - Eta_e(2)) < 100*eps);
                        thetaRange = [theta_x1 theta_x2];
                    case 'West'
                        indices = ~(abs(xi_x - Xi_e(1)) < 100*eps);
                        thetaRange = [theta_x2 theta_x3];
                        indices2 = theta_x3 < 0;
                        thetaRange(indices2,:) = [theta_x2(indices2) theta_x3(indices2)+2*pi];
                end
                for gp = 1:size(W2D_2,1)
                    pt = Q2D_2(gp,:);
                    wt = W2D_2(gp);

                    rho_t = parent2ParametricSpace([0, 1],   pt(1));
%                     theta = parent2ParametricSpace(thetaRange,pt(2));
                    theta = (pt(2)*(thetaRange(:,2) - thetaRange(:,1)) + (thetaRange(:,2) + thetaRange(:,1)))*0.5; % = xiE(1) + (xiTilde+1)*0.5*(xiE(2) - xiE(1)); 
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
                    rho = rho_hat*rho_t;

                    xi_t  = xi_x_t + rho.*cos(theta);
                    eta_t = eta_x_t + rho.*sin(theta);

                    xi = parent2ParametricSpace(Xi_e, xi_t);
                    eta = parent2ParametricSpace(Eta_e, eta_t);

                    [y,u_y,dxdxi,dxdeta] = numericalSolEvalVectorized(xi,eta,p_xi,p_eta,Xi,Eta,wgts,pts,U_sctr);

                    crossProd = cross(dxdxi,dxdeta);
                    J_1 = norm(crossProd);
                    ny = crossProd/J_1;

                    J_3 = rho;
                    J_4 = rho_hat;
                    J_5 = 0.25*(thetaRange(:,2)-thetaRange(:,1));
                    fact = J_1*J_2.*J_3.*J_4.*J_5*wt;
                    u_y = u_y + p_inc(y);
                    xmy = x-y;
                    r = norm2(xmy);
                    dPhi_0dny = Phi_0(r)./r.^2.*              sum(xmy.*ny,2);
                    dPhi_kdny = Phi_k(r)./r.^2.*(1 - 1i*k*r).*sum(xmy.*ny,2);
                    f = f + (dPhi_kdny.*u_y - dPhi_0dny.*u_x).*fact;
                end
            end
        else
            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);
                wt = W2D(gp);

                xi  = parent2ParametricSpace(Xi_e, pt(1));
                eta = parent2ParametricSpace(Eta_e,pt(2));
                
                [y,u_y,dxdxi,dxdeta] = numericalSolEvalVectorized(xi,eta,p_xi,p_eta,Xi,Eta,wgts,pts,U_sctr);

                crossProd = cross(dxdxi,dxdeta);
                J_1 = norm(crossProd);
                ny = crossProd.'/J_1;

                fact = J_1*J_2*wt;

                u_y = u_y + p_inc(y);
                xmy = x-repmat(y,npts^2,1);
                r = norm2(xmy);
                dPhi_0dny = Phi_0(r)./r.^2.*              (xmy*ny);
                dPhi_kdny = Phi_k(r)./r.^2.*(1 - 1i*k*r).*(xmy*ny);
                f = f + (dPhi_kdny*u_y - dPhi_0dny.*u_x)*fact;
            end
        end
    end
    X = reshape(x(:,1),npts,npts);
    Y = reshape(x(:,2),npts,npts);
    Z = reshape(x(:,3),npts,npts);
    C = reshape(f,npts,npts);
%     C = reshape(u_x-(analytic(x)+p_inc(x)),npts,npts);
%     C = reshape(u_x,npts,npts);
%     C = reshape(analytic(x),npts,npts);
%     surf(X,Y,Z, abs(C), 'EdgeColor','none','LineStyle','none')
    surf(X,Y,Z, log10(abs(C)), 'EdgeColor','none','LineStyle','none')
    hold on
    colorbar 
    axis off
%     caxis([-3,0])
    axis equal
    
    set(gca, 'Color', 'none');
%     title(['Fluid 3D NURBS geometry. Mesh ' num2str(M)])
%     view(-70,30)
%     view(120,10)
    view(18,10)
%     view(0,0)
    camproj('perspective')
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    lineColor = 'black';
    lineWidth = 0.5;
    plot3(X(:,1),Y(:,1),Z(:,1),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(:,end),Y(:,end),Z(:,end),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(1,:),Y(1,:),Z(1,:),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(end,:),Y(end,:),Z(end,:),'color',lineColor,'LineWidth',lineWidth)
    drawnow
end
