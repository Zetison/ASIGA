function maxU = createParaviewFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options, e3Dss_options, rho)

maxU = NaN;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
Eps = 1e4*eps;

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
noDofs = varCol.noDofs;
patches = varCol.patches;
d = varCol.dimension;
if isfield('varCol','omega')
    omega = varCol.omega;
else
    omega = NaN;
end
isOuterDomain = varCol.isOuterDomain;
if isOuterDomain
    model = varCol.model;
end

if d == 1
    if isfield(varCol, 'analytic')
        analytic = varCol.analytic;
    else
        analytic = 0;
    end
    if isfield(varCol, 'gAnalytic')
        gAnalytic = varCol.gAnalytic;
    else
        gAnalytic = 0;
    end
    
    if isOuterDomain && ~strcmp(model,'PS')
        p_inc = varCol.p_inc;
        gp_inc = varCol.gp_inc;
    end
end
stringShift = 40;
if d == 3
    C = varCol.C;
    Ux = U(1:d:noDofs);
    Uy = U(2:d:noDofs);
    Uz = U(3:d:noDofs);
    U = [Ux, Uy, Uz];
else
%     varCol.dofsToRemove = varCol.dofsToRemove_old;
%     tic
%     fprintf(['\n%-' num2str(stringShift) 's'], '    Performing least squares ... ')
%     dU = leastSquares(varCol,U,'gradient');
%     dU = leastSquares(varCol,U,'scalar');
%     fprintf('using %12f seconds.', toc)
    C = 0;
end
type = patches{1}.nurbs.type;
switch type
    case '3Dsurface'
        container = cell(1,noElems);
        parfor e = 1:noElems
            patch = pIndex(e);
            Xi = knotVecs{patch}{1};
            Eta = knotVecs{patch}{2};

            idXi = index(e,1);
            idEta = index(e,2);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New
            Usctr = U(sctr,:);
            
            noXiKnots = 2+extraXiPts;
            noEtaKnots = 2+extraEtaPts;
            noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
            visElements_e = zeros(noVisElems,4);
            eVis = 1;
            for j = 1:noEtaKnots-1
                for i = 1:noXiKnots-1
                    visElements_e(eVis,1) = i   +   (j-1)*noXiKnots;
                    visElements_e(eVis,2) = i+1 +   (j-1)*noXiKnots;
                    visElements_e(eVis,3) = i+1 +       j*noXiKnots;
                    visElements_e(eVis,4) = i   +       j*noXiKnots;
                    eVis = eVis + 1;
                end
            end
            noNodes = noXiKnots*noEtaKnots;
            scalarField_e = zeros(noNodes, 1);    
            testField_e = zeros(noNodes, 1);            
            xiKnots = linspace(Xi_e(1)+Eps,Xi_e(2)-Eps,noXiKnots);        
            etaKnots = linspace(Eta_e(1)+Eps,Eta_e(2)-Eps,noEtaKnots);
            counter = 1;
            nodes_e = zeros(noNodes,3);
            for j = 1:noEtaKnots
                eta = etaKnots(j);
                for i = 1:noXiKnots
                    xi = xiKnots(i);
                    [u, v] = numericalSolEval_final_surf(xi, eta, p_xi, p_eta, Xi, Eta, wgts, pts, Usctr);
                    R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
                    scalarField_e(counter) = u;
%                     testField_e(counter) = sum(R.^2);
                    testField_e(counter) = sqrt(sum(wgts.^2)/numel(wgts));
%                     testField_e(counter) = sum(R);
                    nodes_e(counter,:) = v;
                    counter = counter + 1;
                end
            end
            container{e}.nodes = nodes_e;
            container{e}.noNodes = noNodes;
            container{e}.noVisElems = noVisElems;
            container{e}.scalarField = scalarField_e;
            container{e}.testField = testField_e;
            container{e}.visElements = visElements_e;
        end
        noNodes = 0;
        noVisElems = 0;
        for e = 1:noElems
            noNodes = noNodes + container{e}.noNodes;
            noVisElems = noVisElems + container{e}.noVisElems;
        end
        visElements = zeros(noVisElems,4);
        scalarField = zeros(noNodes,1);
        testField = zeros(noNodes,1);
        nodes = zeros(noNodes,3);
        nodesCount = 0;
        count_vis = 0;
        for e = 1:noElems
            visElements(count_vis+1:count_vis+container{e}.noVisElems,:) = nodesCount + container{e}.visElements;
            nodes(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.nodes;
            scalarField(nodesCount+1:nodesCount+container{e}.noNodes) = container{e}.scalarField;
            testField(nodesCount+1:nodesCount+container{e}.noNodes) = container{e}.testField;
            nodesCount = nodesCount + container{e}.noNodes;
            count_vis = count_vis + container{e}.noVisElems;
        end
    case '3Dvolume'
        p_zeta = varCol.degree(3); % assume p_eta is equal in all patches
        elRangeZeta = varCol.elRange{3};
        n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);
        noXiKnots = 2+extraXiPts;
        noEtaKnots = 2+extraEtaPts;
        noZetaKnots = 2+extraZetaPts;
        noNodes = noXiKnots*noEtaKnots*noZetaKnots;   
        noVisElems  = (noXiKnots-1)*(noEtaKnots-1)*(noZetaKnots-1);
        container = cell(1,noElems);
        for e = 1:noElems
            container{e}.nodes = zeros(noNodes,3);
            container{e}.noNodes = noNodes;
            container{e}.noVisElems = noVisElems;
            container{e}.displacement = zeros(noNodes,3);
            container{e}.density = zeros(noNodes,1);
            container{e}.visElements = zeros(noVisElems,8);
            container{e}.strain = zeros(noNodes,6);
        end
%         for e = 1:noElems
        parfor e = 1:noElems
            patch = pIndex(e);
            Xi = knotVecs{patch}{1};
            Eta = knotVecs{patch}{2};
            Zeta = knotVecs{patch}{3};

            idXi = index(e,1);
            idEta = index(e,2);
            idZeta = index(e,3);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);
            Zeta_e = elRangeZeta(idZeta,:);

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New
            Usctr = U(sctr,:);
            rhosctr = rho(sctr);
            
            visElements_e = zeros(noVisElems,8);
            eVis = 1;
            for k = 1:noZetaKnots-1
                for j = 1:noEtaKnots-1
                    for i = 1:noXiKnots-1
                        visElements_e(eVis,1) = i   +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                        visElements_e(eVis,2) = i+1 +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                        visElements_e(eVis,3) = i+1 +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                        visElements_e(eVis,4) = i   +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                        visElements_e(eVis,5) = i   +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
                        visElements_e(eVis,6) = i+1 +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
                        visElements_e(eVis,7) = i+1 +       j*noXiKnots +       k*noEtaKnots*noXiKnots;
                        visElements_e(eVis,8) = i   +       j*noXiKnots +       k*noEtaKnots*noXiKnots;

                        eVis = eVis + 1;
                    end
                end
            end  
            xiKnots = linspace(Xi_e(1)+Eps,Xi_e(2)-Eps,noXiKnots);        
            etaKnots = linspace(Eta_e(1)+Eps,Eta_e(2)-Eps,noEtaKnots);
            zetaKnots = linspace(Zeta_e(1)+Eps,Zeta_e(2)-Eps,noZetaKnots);
            [xi,eta,zeta] = ndgrid(xiKnots,etaKnots,zetaKnots);
            xi = xi(:);
            eta = eta(:);
            zeta = zeta(:);
            
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
            
            container{e}.nodes = R*pts;
            container{e}.displacement = R*Usctr;
            container{e}.density = R*rhosctr;
            container{e}.visElements = visElements_e;
            
            indices = abs(J_1) < 100*Eps;
            if any(indices)
                digits(100)
                [~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVec(vpa(xi(indices)), vpa(eta(indices)), vpa(zeta(indices)), p_xi, p_eta, p_zeta, vpa(Xi), vpa(Eta), vpa(Zeta), vpa(wgts));
                dXdxi = dRdxi*vpa(pts);
                dXdeta = dRdeta*vpa(pts);
                dXdzeta = dRdzeta*vpa(pts);
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
                dRdx(indices,:) = repmat(Jinv1(:,1),1,n_en).*dRdxi + repmat(Jinv1(:,2),1,n_en).*dRdeta + repmat(Jinv1(:,3),1,n_en).*dRdzeta;
                dRdy(indices,:) = repmat(Jinv2(:,1),1,n_en).*dRdxi + repmat(Jinv2(:,2),1,n_en).*dRdeta + repmat(Jinv2(:,3),1,n_en).*dRdzeta;
                dRdz(indices,:) = repmat(Jinv3(:,1),1,n_en).*dRdxi + repmat(Jinv3(:,2),1,n_en).*dRdeta + repmat(Jinv3(:,3),1,n_en).*dRdzeta;
            end
            dudx = dRdx*Usctr;
            dudy = dRdy*Usctr;
            dudz = dRdz*Usctr;
            container{e}.strain = calculateStrainVectorVec(dudx, dudy, dudz);
        end
        noNodes = 0;
        noVisElems = 0;
        for e = 1:noElems
            noNodes = noNodes + container{e}.noNodes;
            noVisElems = noVisElems + container{e}.noVisElems;
        end
        visElements = zeros(noVisElems,8);
        displacement = zeros(noNodes,3);
        density = zeros(noNodes,1);
        strain = zeros(noNodes,6);
        nodes = zeros(noNodes,3);
        nodesCount = 0;
        count_vis = 0;
        for e = 1:noElems
            visElements(count_vis+1:count_vis+container{e}.noVisElems,:) = nodesCount + container{e}.visElements;
            nodes(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.nodes;
            displacement(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.displacement;
            density(nodesCount+1:nodesCount+container{e}.noNodes) = container{e}.density;
            strain(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.strain;
            nodesCount = nodesCount + container{e}.noNodes;
            count_vis = count_vis + container{e}.noVisElems;
        end
end
data.nodes = nodes;
data.visElements = visElements;
data.omega = omega;

% tic
% fprintf(['\n%-' num2str(stringShift) 's'], '    Computing error/storing data ... ')
if isOuterDomain
    data.P_inc = real(makeDynamic(p_inc(nodes), options, omega)); 
    if strcmp(varCol.method,'BEM')
        totField = scalarField;
        if isOuterDomain && ~strcmp(model,'PS')
            scalarField = scalarField - data.P_inc;
        end
        data.errorFunc = errorFunc;
    else   
        totField = scalarField + p_inc(nodes);
    end
    if ~isempty(e3Dss_options)
        
        if strcmp(varCol.applyLoad, 'radialPulsation')
            data_e3Dss.p = varCol.analytic(nodes);
            dp = varCol.gAnalytic(nodes);
            data_e3Dss.dpdx = dp(:,1);
            data_e3Dss.dpdy = dp(:,2);
            data_e3Dss.dpdz = dp(:,3);
            data.Error = abs(data_e3Dss.p-scalarField)./abs(data_e3Dss.p);
        else
            data_e3Dss = e3Dss(nodes,e3Dss_options);
        end
        analytic_v = data_e3Dss(1).p;
        if strcmp(type, '3Dvolume')
            p = data_e3Dss(1).p;
            dp = [data_e3Dss(1).dpdx, data_e3Dss(1).dpdy, data_e3Dss(1).dpdz];
            p_e = p-scalarField;
%             dp_e = dp-gScalarField_p;

            p2 = abs(p).^2;
            dp2 = sum(abs(dp).^2,2);
            p_e2 = abs(p_e).^2;
            dp_e2 = sum(abs(dp_e).^2,2);
            
            k = varCol.k;
            data.Error = sqrt(p_e2/max(p2));
%             data.ErrorGrad = sqrt(dp_e2/max(dp2));
%             data.ErrorEnergy = sqrt((dp_e2 + k^2*p_e2)/max(dp2 + k^2*p2));
        end
        data.analytic = real(makeDynamic(analytic_v, options, omega)); 
    end
    data.testField = testField;
%     data.testFun = varCol.testFun(nodes); 
    data.scalarField = real(makeDynamic(scalarField, options, omega)); 
    data.totField = real(makeDynamic(totField, options, omega));
    data.totFieldAbs = abs(makeDynamic(totField, options, omega));
%     data.testFun = abs(norm2(gScalarField2-gScalarField)/max(norm2(gScalarField)));
%     data.testFun = abs(norm2(scalarField-gScalarField)/max(abs(norm2(scalarField))));
%     data.testFun = abs((norm2(gScalarField2-gAnalytic_v,0) - k(1)^2*abs(scalarField-analytic_v))/max(norm2(gAnalytic_v,0) - k(1)^2*abs(analytic_v)));
    
    if strcmp(type, '3Dvolume')
        if d == 1 && isOuterDomain
            rho_f = varCol.rho_f;
            displacement = (gScalarField_p+gp_inc(nodes))/(rho_f*omega^2);
        end
    end    
elseif d == 1
    data.totField = real(makeDynamic(scalarField, options, omega));
    data.totFieldAbs = abs(makeDynamic(scalarField, options, omega));
    rho_f = varCol.rho_f; 
    displacement = gScalarField_p/(rho_f*omega^2);  
    if ~isempty(e3Dss_options)
        data_e3Dss = e3Dss({[],[],nodes},e3Dss_options);
        p = data_e3Dss(2).p;
        dp = [data_e3Dss(2).dpdx, data_e3Dss(2).dpdy, data_e3Dss(2).dpdz];
        p_e = p-scalarField;
        dp_e = dp-gScalarField_p;

        p2 = abs(p).^2;
        dp2 = sum(abs(dp).^2,2);
        p_e2 = abs(p_e).^2;
        dp_e2 = sum(abs(dp_e).^2,2);

        k = varCol.k;
        data.Error = sqrt(p_e2/max(p2));
        data.ErrorGrad = sqrt(dp_e2/max(dp2));
        data.ErrorEnergy = sqrt((dp_e2 + k^2*p_e2)/max(dp2 + k^2*p2));
    end
else
    if ~isempty(e3Dss_options)
        e3Dss_options.calc_sigma_xx = 1;
        e3Dss_options.calc_sigma_yy = 1; 
        e3Dss_options.calc_sigma_zz = 1;
        e3Dss_options.calc_sigma_yz = 1;
        e3Dss_options.calc_sigma_xz = 1;
        e3Dss_options.calc_sigma_xy = 1; 
        e3Dss_options.calc_u_x = 1; 
        e3Dss_options.calc_u_y = 1; 
        e3Dss_options.calc_u_z = 1; 
        e3Dss_options.calc_du_xdx = 1; 
        e3Dss_options.calc_du_xdy = 1; 
        e3Dss_options.calc_du_xdz = 1; 
        e3Dss_options.calc_du_ydx = 1; 
        e3Dss_options.calc_du_ydy = 1; 
        e3Dss_options.calc_du_ydz = 1; 
        e3Dss_options.calc_du_zdx = 1; 
        e3Dss_options.calc_du_zdy = 1; 
        e3Dss_options.calc_du_zdz = 1; 
        data_e3Dss = e3Dss({[], nodes},e3Dss_options);

        u = [data_e3Dss(1).u_x,data_e3Dss(1).u_y,data_e3Dss(1).u_z];                
        u_e = u-displacement;
                
        u2 = sum(abs(u).^2,2);
        u_e2 = sum(abs(u_e).^2,2);
        
        
        
        sigma = [data_e3Dss(1).sigma_xx,data_e3Dss(1).sigma_yy,data_e3Dss(1).sigma_zz,data_e3Dss(1).sigma_yz,data_e3Dss(1).sigma_xz,data_e3Dss(1).sigma_xy];
        C = varCol.C;
        
        strain_vec = (C\sigma.').';
        
        strain_e = strain_vec-strain_h;
        
        eCe = real(sum((strain_e*C).*conj(strain_e),2)); % the usage of real() is to remove machine epsilon imaginary part
        uCu = real(sum((strain_vec*C).*conj(strain_vec),2)); % the usage of real() is to remove machine epsilon imaginary part
        
        rho_s = varCol.rho_s; 
        data.Error = sqrt(u_e2/max(u2));
        data.ErrorGrad = sqrt(eCe/max(uCu));
        data.ErrorEnergy = sqrt((eCe + rho_s*omega^2*u_e2)/max(uCu + rho_s*omega^2*u2));
    end
    data.stress = real(makeDynamic(strain, options, omega));
    data.displacement = real(makeDynamic(displacement, options, omega));
    data.density = real(makeDynamic(density, options, omega));
end
% fprintf('using %12f seconds.', toc)
% data.displacement = real(makeDynamic(displacement, options, omega));
% keyboard
% tic
% fprintf(['\n%-' num2str(stringShift) 's'], '    Creating VTK-file ... ')
% keyboard

makeVTKfile(data, options);
% fprintf('using %12f seconds.', toc)
