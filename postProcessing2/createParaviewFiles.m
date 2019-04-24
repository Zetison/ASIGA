function maxU = createParaviewFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options, plotMesh, e3Dss_options)

maxU = NaN;

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
noDofs = varCol.noDofs;
patches = varCol.patches;
d = varCol.dimension;
omega = varCol.omega;
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
            
%             nurbs = patches.nurbs;
%             visObj_p = buildVisualizationMesh(nurbs, [extraXiPts, extraEtaPts]);
%             nodes_p = visObj_p.nodes;
%             noNodes = visObj_p.noNodes;
%             visElements = visObj_p.visElements;
%             noXiKnots = visObj_p.noXiKnots;
%             noEtaKnots = visObj_p.noEtaKnots;
%             XiVec = visObj_p.XiVec;
%             EtaVec = visObj_p.EtaVec;
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
            xiKnots = linspace(Xi_e(1),Xi_e(2)-eps,noXiKnots);        
            etaKnots = linspace(Eta_e(1),Eta_e(2)-eps,noEtaKnots);
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
        tic
        fprintf(['\n%-' num2str(stringShift) 's'], '    Creating post-processing mesh ... ')
        Zeta = varCol.nurbs.knots{3};
        p_zeta = varCol.nurbs.degree(3);
        
        visObj_p = buildVisualizationMesh(nurbs, [extraXiPts, extraEtaPts, extraZetaPts]);
        fprintf('using %12f seconds.', toc)
        
        
%         keyboard
        nodes = visObj_p.nodes;
        visElements = visObj_p.visElements;

        tic
        fprintf(['\n%-' num2str(stringShift) 's'], '    Compute solution ... ')
        multiValuedNodes = 0;
        if multiValuedNodes
            nodes = visObj_p.nodes;
            elRangeXi = visObj_p.elRangeXi;
            elRangeEta = visObj_p.elRangeEta;
            elRangeZeta = visObj_p.elRangeZeta;
            index = visObj_p.index;
            d_test = 3;
    %         d_test = 1;
            noVisElems = size(visElements,1);
            
            scalarField = zeros(8, noVisElems);
%             gScalarField = zeros(8, noVisElems, d_test);
            gScalarField_p = zeros(8, noVisElems, 3);
            strain_h = zeros(8, noVisElems, 6);
            displacement = zeros(8, noVisElems, 3);
            parfor e = 1:noVisElems
%             for e = 1:noVisElems
                idXi = index(e,1);
                idEta = index(e,2);
                idZeta = index(e,3);
                
                scalarField_temp = zeros(8,1);
%                 gScalarField_temp2 = zeros(8, d_test);
                gScalarField_temp = zeros(8, 3);
                displacement_temp = zeros(8, 3);
                strain_temp = zeros(8, 6);
                Xi_e = elRangeXi(idXi,:);
                Eta_e = elRangeEta(idEta,:);
                Zeta_e = elRangeZeta(idZeta,:);
                rot = [1, 4; 
                       2, 3];
                for l = 1:2
                    zeta = Zeta_e(l) + eps*(3-2*l);
                    for j = 1:2
                        eta = Eta_e(j) + eps*(3-2*j);
                        for i = 1:2
                            xi = Xi_e(i) + eps*(3-2*i);
                            count = 4*(l-1) + rot(i,j);
%                             count = 4*(l-1) + 2*(j-1) + i;
                            if d == 1
                %                 u = numericalSolEval_final(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights, controlPts, U);
%                                 u = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
%                                 du = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, dU);
                                
%                                 gScalarField_temp(count,:) = du; %[dudx, dudy, dudz];
                                [u, ~, dudx, dudy, dudz] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                                scalarField_temp(count) = u;
                                gScalarField_temp(count,:) = [dudx, dudy, dudz];
                %                 keyboard
                            elseif d == 3
                                [u, ~, dudx, dudy, dudz] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                                displacement_temp(count,:) = u;

                                strain_temp(count,:) = calculateStrainVector(dudx, dudy, dudz);
                            end
                        end
                    end
                end
                scalarField(:,e) = scalarField_temp;
%                 gScalarField(:,e,:) = gScalarField_temp2;
                gScalarField_p(:,e,:) = gScalarField_temp;
                displacement(:,e,:) = displacement_temp;
                strain_h(:,e,:) = strain_temp;
                du_h(:,e,:) = du_h_temp;
            end
            temp = visElements.';
            nodes = nodes(temp(:),:);
            visElements = reshape((1:numel(visElements)).',8,[]).';
            gScalarField_p = reshape(gScalarField_p, noVisElems*8, d_test);
            gScalarField_p = reshape(gScalarField_p, noVisElems*8, 3);
            scalarField = reshape(scalarField, noVisElems*8, 1);
            % errorFunc = reshape(errorFunc, noVisElems*8, 1);
            displacement = reshape(displacement, noVisElems*8, 3);
            strain_h = reshape(strain_h, noVisElems*8, 6);
        else
            noNodes = visObj_p.noNodes;
            noXiKnots = visObj_p.noXiKnots;
            noEtaKnots = visObj_p.noEtaKnots;
            noZetaKnots = visObj_p.noZetaKnots;
            XiVec = visObj_p.XiVec;
            EtaVec = visObj_p.EtaVec;
            ZetaVec = visObj_p.ZetaVec;
            scalarField = zeros(noEtaKnots*noZetaKnots, noXiKnots);
            
            gScalarField_p = zeros(noEtaKnots*noZetaKnots, noXiKnots, 3);
            strain_h = zeros(noEtaKnots*noZetaKnots, noXiKnots, 6);
            displacement = zeros(noEtaKnots*noZetaKnots, noXiKnots, 3);
            
            parfor i = 1:noXiKnots
%             for i = 1:noXiKnots
                scalarField_temp = zeros(noEtaKnots*noZetaKnots, 1);
%                 gScalarField_temp2 = zeros(noEtaKnots*noZetaKnots, d_test);
                gScalarField_temp = zeros(noEtaKnots*noZetaKnots, 3);
                displacement_temp = zeros(noEtaKnots*noZetaKnots, 3);
                strain_temp = zeros(noEtaKnots*noZetaKnots, 6);
                count = 1;
                xi = XiVec(i);
                for j = 1:noEtaKnots
                    eta = EtaVec(j);
                    for l = 1:noZetaKnots
                        zeta = ZetaVec(l);    


                        if d == 1
            %                 u = numericalSolEval_final(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights, controlPts, U);
%                             u = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
%                             du = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, dU);
%                             gScalarField_temp2(count,:) = du; %[dudx, dudy, dudz];
                            [u, ~, dudx, dudy, dudz] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                            scalarField_temp(count) = u;
                            gScalarField_temp(count,:) = [dudx, dudy, dudz];
            %                 keyboard
                        elseif d == 3
                            [u, ~, dudx, dudy, dudz] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                            displacement_temp(count,:) = u;

                            strain_temp(count,:) = calculateStrainVector(dudx, dudy, dudz);
                        end
                        count = count + 1;
                    end
                end
                scalarField(:,i) = scalarField_temp;
%                 gScalarField(:,i,:) = gScalarField_temp2;
                gScalarField_p(:,i,:) = gScalarField_temp;
                displacement(:,i,:) = displacement_temp;
                strain_h(:,i,:) = strain_temp;
            end
            gScalarField_p = reshape(gScalarField_p, noZetaKnots, noEtaKnots, noXiKnots, 3);
            gScalarField_p = permute(gScalarField_p, [3,2,1,4]);
            gScalarField_p = reshape(gScalarField_p, noNodes, 3);
            
            scalarField = reshape(scalarField, noZetaKnots, noEtaKnots, noXiKnots);
            scalarField = permute(scalarField, [3,2,1]);
            scalarField = reshape(scalarField, noNodes, 1);
            
            displacement = reshape(displacement, noZetaKnots, noEtaKnots, noXiKnots, 3);
            displacement = permute(displacement, [3,2,1,4]);
            displacement = reshape(displacement, noNodes, 3);
            
            strain_h = reshape(strain_h, noZetaKnots, noEtaKnots, noXiKnots, 6);
            strain_h = permute(strain_h, [3,2,1,4]);
            strain_h = reshape(strain_h, noNodes, 6);
        end
        fprintf('using %12f seconds.', toc)
end
% data.gradient = real(makeDynamic(gScalarField_p, options, omega)); 
data.nodes = nodes;
data.visElements = visElements;
data.omega = omega;

tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Computing error/storing data ... ')
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
    data.stress = real(makeDynamic(strain_h, options, omega));
end
fprintf('using %12f seconds.', toc)
% data.displacement = real(makeDynamic(displacement, options, omega));
% keyboard
tic
fprintf(['\n%-' num2str(stringShift) 's'], '    Creating VTK-file ... ')
% keyboard

makeVTKfile(data, options);
fprintf('using %12f seconds.', toc)
