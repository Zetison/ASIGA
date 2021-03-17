function maxU = createParaviewFiles(varargin)
varCol = varargin{1};
options = struct('para_options', 0,...
                 'rho',NaN);
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
analyticSolutionExist = varCol{1}.analyticSolutionExist;
splitExteriorFields = varCol{1}.splitExteriorFields;
noDomains = numel(varCol);
nodes = cell(1,noDomains);
visElements = cell(1,noDomains);
scalarField = cell(1,noDomains);
para = cell(1,noDomains);
displacement = cell(1,noDomains);
density = cell(1,noDomains);
strain_h = cell(1,noDomains);
gScalarField_p = cell(1,noDomains);
stringShift = 40;
fprintf(['\n%-' num2str(stringShift) 's'], '    Evaluating solution ... ')
tic
for i_v = 1:noDomains
    para{i_v} = options.para_options;
    if varCol{i_v}.boundaryMethod
        celltype = 'VTK_QUAD';
    else
        celltype = 'VTK_HEXAHEDRON';
    end
    isOuterDomain = i_v == 1;
    varCol{i_v}.isOuterDomain = isOuterDomain;
    d = varCol{i_v}.dimension;
    isSolid = d == 3;

    para{i_v}.celltype = celltype;
    para{i_v}.plotP_inc = options.para_options.plotP_inc && ~isSolid && splitExteriorFields && isOuterDomain;
    para{i_v}.plotScalarField = options.para_options.plotScalarField && ~isSolid && isOuterDomain;
    para{i_v}.plotScalarFieldAbs = options.para_options.plotScalarField && ~isSolid && isOuterDomain;
    para{i_v}.plotTotField = options.para_options.plotTotField && ~isSolid && splitExteriorFields; 
    para{i_v}.plotTotFieldAbs = options.para_options.plotTotFieldAbs && ~isSolid && splitExteriorFields; 
    para{i_v}.plotAnalytic = options.para_options.plotAnalytic && analyticSolutionExist; 
    para{i_v}.computeGrad = options.para_options.computeGrad && ~(varCol{i_v}.boundaryMethod && ~isSolid);
    para{i_v}.plotError = options.para_options.plotError && (analyticSolutionExist && ~options.para_options.plotTimeOscillation); 
    para{i_v}.plotErrorGrad = para{i_v}.plotError; 
    para{i_v}.plotErrorEnergy = para{i_v}.plotError; 
    
    para{i_v}.plotVonMisesStress = options.para_options.plotVonMisesStress && isSolid;
    para{i_v}.plotStressXX = options.para_options.plotStressXX && isSolid;
    para{i_v}.plotStressYY = options.para_options.plotStressYY && isSolid;
    para{i_v}.plotStressZZ = options.para_options.plotStressZZ && isSolid;
    para{i_v}.plotStressYZ = options.para_options.plotStressYZ && isSolid;
    para{i_v}.plotStressXZ = options.para_options.plotStressXZ && isSolid;
    para{i_v}.plotStressXY = options.para_options.plotStressXY && isSolid;

    U = varCol{i_v}.U;
    rho = options.rho;
    if isnan(options.rho)
        rho = zeros(size(U,1),1);
        withDensity = false;
    else
        withDensity = true;

    end
    maxU = NaN;

    degree = varCol{i_v}.degree;

    extraXiPts = options.para_options.extraXiPts;
    extraEtaPts = options.para_options.extraEtaPts;
    extraZetaPts = options.para_options.extraZetaPts;

    index = varCol{i_v}.index;
    noElems = varCol{i_v}.noElems;
    elRange = varCol{i_v}.elRange;
    element = varCol{i_v}.element;
    element2 = varCol{i_v}.element2;
    weights = varCol{i_v}.weights;
    controlPts = varCol{i_v}.controlPts;
    knotVecs = varCol{i_v}.knotVecs;
    pIndex = varCol{i_v}.pIndex;
    noDofs = varCol{i_v}.noDofs;
    patches = varCol{i_v}.patches;
    d = varCol{i_v}.dimension;
    if isfield(varCol{i_v},'omega')
        omega = varCol{i_v}.omega;
    else
        omega = NaN;
    end

    if d == 3
        C = varCol{i_v}.C;
        Ux = U(1:d:noDofs);
        Uy = U(2:d:noDofs);
        Uz = U(3:d:noDofs);
        U = [Ux, Uy, Uz];
    else
    %     varCol{i_v}.dofsToRemove = varCol{i_v}.dofsToRemove_old;
    %     tic
    %     fprintf(['\n%-' num2str(stringShift) 's'], '    Performing least squares ... ')
    %     dU = leastSquares(varCol{i_v},U,'gradient');
    %     dU = leastSquares(varCol{i_v},U,'scalar');
    %     fprintf('using %12f seconds.', toc)
        C = 0;
    end
    n_en = prod(degree+1);
    d_p = patches{1}.nurbs.d_p;
    if options.para_options.plotDisplacementVectors && strcmp(varCol{1}.BC,'SHBC') && varCol{1}.boundaryMethod
        warning('Displacement (gradient of pressure) may not be plotted for SHBC and isBoundaryMethod')
    end
    switch d_p
        case 2
            noXiKnots = 2+extraXiPts;
            noEtaKnots = 2+extraEtaPts;
            noNodes = noXiKnots*noEtaKnots;   
            noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
            container = cell(1,noElems);
            for e = 1:noElems
                container{e}.nodes = zeros(noNodes,3);
                container{e}.noNodes = noNodes;
                container{e}.noVisElems = noVisElems;
                container{e}.displacement{i_v} = zeros(noNodes,3);
                container{e}.scalarField{i_v} = zeros(noNodes,3);
                container{e}.density{i_v} = zeros(noNodes,1);
                container{e}.visElements = zeros(noVisElems,8);
                container{e}.strain = zeros(noNodes,6);
            end
%             for e = 1:noElems
            parfor e = 1:noElems
                patch = pIndex(e);
                knots = knotVecs{patch};

                Xi_e = zeros(d_p,2);
                for i = 1:d_p
                    Xi_e(i,:) = elRange{i}(index(e,i),:);
                end

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
                xiKnots = linspace(Xi_e(1,1),Xi_e(1,2)-eps,noXiKnots);        
                etaKnots = linspace(Xi_e(2,1),Xi_e(2,2)-eps,noEtaKnots);   
                [xi,eta] = ndgrid(xiKnots,etaKnots);
                xi = xi(:);
                eta = eta(:);

                I = findKnotSpans(degree, [xi(1),eta(1)], knots);
                R = NURBSbasis(I,[xi, eta], degree, knots, wgts);
    %             dXdxi = R{2}*pts;
    %             dXdeta = R{3}*pts;

                nodes_e = R{1}*pts;
                scalarField_e = R{1}*Usctr;

                container{e}.nodes = nodes_e;
                container{e}.noNodes = noNodes;
                container{e}.noVisElems = noVisElems;
                container{e}.scalarField{i_v} = scalarField_e;
                container{e}.visElements = visElements_e;
            end
            noNodes = 0;
            noVisElems = 0;
            for e = 1:noElems
                noNodes = noNodes + container{e}.noNodes;
                noVisElems = noVisElems + container{e}.noVisElems;
            end
            visElements{i_v} = zeros(noVisElems,4);
            scalarField{i_v} = zeros(noNodes,d);
            nodes{i_v} = zeros(noNodes,3);
            nodesCount = 0;
            count_vis = 0;
            for e = 1:noElems
                visElements{i_v}(count_vis+1:count_vis+container{e}.noVisElems,:) = nodesCount + container{e}.visElements;
                nodes{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.nodes;
                scalarField{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.scalarField{i_v};
                displacement{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.displacement{i_v};
                nodesCount = nodesCount + container{e}.noNodes;
                count_vis = count_vis + container{e}.noVisElems;
            end
        case 3
            Eps = 1e-6;
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
                container{e}.displacement{i_v} = zeros(noNodes,3);
                container{e}.scalarField{i_v} = zeros(noNodes,1);
                container{e}.gScalarField_p{i_v} = zeros(noNodes,3);
                container{e}.density{i_v} = zeros(noNodes,1);
                container{e}.visElements = zeros(noVisElems,8);
                container{e}.strain = zeros(noNodes,6);
            end
%             for e = 1:noElems
            parfor e = 1:noElems
                patch = pIndex(e);
                knots = knotVecs{patch};

                Xi_e = zeros(d_p,2);
                for i = 1:d_p
                    Xi_e(i,:) = elRange{i}(index(e,i),:);
                end

                sctr = element(e,:);
                pts = controlPts(sctr,:);
                wgts = weights(element2(e,:),:); % New
                Usctr = U(sctr,:);

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
                xiKnots = linspace(Xi_e(1,1),Xi_e(1,2)-eps,noXiKnots);        
                etaKnots = linspace(Xi_e(2,1),Xi_e(2,2)-eps,noEtaKnots);    
                zetaKnots = linspace(Xi_e(3,1),Xi_e(3,2)-eps,noZetaKnots);
                [xi,eta,zeta] = ndgrid(xiKnots,etaKnots,zetaKnots);
                xi = xi(:);
                eta = eta(:);
                zeta = zeta(:);

                I = findKnotSpans(degree, [xi(1),eta(1),zeta(1)], knots);
                R = NURBSbasis(I,[xi, eta, zeta], degree, knots, wgts);
                dXdxi = R{2}*pts;
                dXdeta = R{3}*pts;
                dXdzeta = R{4}*pts;
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
                dRdx = repmat(Jinv1(:,1),1,n_en).*R{2} + repmat(Jinv1(:,2),1,n_en).*R{3} + repmat(Jinv1(:,3),1,n_en).*R{4};
                dRdy = repmat(Jinv2(:,1),1,n_en).*R{2} + repmat(Jinv2(:,2),1,n_en).*R{3} + repmat(Jinv2(:,3),1,n_en).*R{4};
                dRdz = repmat(Jinv3(:,1),1,n_en).*R{2} + repmat(Jinv3(:,2),1,n_en).*R{3} + repmat(Jinv3(:,3),1,n_en).*R{4};

                container{e}.nodes = R{1}*pts;
                container{e}.visElements = visElements_e;
                dudx = dRdx*Usctr;
                dudy = dRdy*Usctr;
                dudz = dRdz*Usctr;
                singularJ = abs(J_1) < Eps;
                dudx(or(singularJ, isnan(dudx))) = 0;
                dudy(or(singularJ, isnan(dudy))) = 0;
                dudz(or(singularJ, isnan(dudz))) = 0;
                if d == 3
                    container{e}.strain = calculateStrainVectorVec(dudx, dudy, dudz);
                    container{e}.displacement{i_v} = R{1}*Usctr;
                else
                    container{e}.scalarField{i_v} = R{1}*Usctr;
                    container{e}.gScalarField_p{i_v} = [dudx,dudy,dudz];
                end
                if withDensity
                    container{e}.density{i_v} = R{1}*rho(sctr);
                end
            end
            noNodes = 0;
            noVisElems = 0;
            for e = 1:noElems
                noNodes = noNodes + container{e}.noNodes;
                noVisElems = noVisElems + container{e}.noVisElems;
            end
            visElements{i_v} = zeros(noVisElems,8);
            displacement{i_v} = zeros(noNodes,3);
            density{i_v} = zeros(noNodes,1);
            strain_h{i_v} = zeros(noNodes,6);
            nodes{i_v} = zeros(noNodes,3);
            gScalarField_p{i_v} = zeros(noNodes,3);
            nodesCount = 0;
            count_vis = 0;
            for e = 1:noElems
                visElements{i_v}(count_vis+1:count_vis+container{e}.noVisElems,:) = nodesCount + container{e}.visElements;
                nodes{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.nodes;
                if d == 3
                    displacement{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.displacement{i_v};
                    strain_h{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.strain;
                else
                    scalarField{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.scalarField{i_v};
                    gScalarField_p{i_v}(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.gScalarField_p{i_v};
                end
                if withDensity
                    density{i_v}(nodesCount+1:nodesCount+container{e}.noNodes) = container{e}.density{i_v};
                end
                nodesCount = nodesCount + container{e}.noNodes;
                count_vis = count_vis + container{e}.noVisElems;
            end
    end
end
fprintf('using %12f seconds.', toc)
if analyticSolutionExist
    layer = varCol{1}.analyticFunctions(nodes);
end

for i_v = 1:numel(varCol)
    fprintf(['\n%-' num2str(stringShift) 's'], '    Computing error/storing data ... ')
    tic
    isOuterDomain = i_v == 1;
    data.nodes = nodes{i_v};
    data.visElements = visElements{i_v};
    data.omega = omega;
    switch varCol{i_v}.media
        case 'fluid'
            if isOuterDomain
                if para{i_v}.plotP_inc
                    data.P_inc = real(makeDynamic(varCol{i_v}.p_inc_(nodes{i_v}), para{i_v}, omega)); 
                end
                if varCol{i_v}.solveForPtot
                    totField = scalarField{i_v};
                    if isOuterDomain && splitExteriorFields
                        scalarField{i_v} = scalarField{i_v} - varCol{i_v}.p_inc_(nodes{i_v});
                    end
                else   
                    if splitExteriorFields
                        totField = scalarField{i_v} + varCol{i_v}.p_inc_(nodes{i_v});
                    else
                        totField = scalarField{i_v};
                    end
                end
                data.scalarField = real(makeDynamic(scalarField{i_v}, para{i_v}, omega)); 
                data.scalarFieldAbs = abs(makeDynamic(scalarField{i_v}, para{i_v}, omega)); 
                if d_p == 3
                    rho_f = varCol{i_v}.rho;
                    if splitExteriorFields
                        gp_inc = [varCol{i_v}.dp_incdx_(nodes{i_v}),varCol{i_v}.dp_incdy_(nodes{i_v}),varCol{i_v}.dp_incdz_(nodes{i_v})];
                        displacement{i_v} = (gScalarField_p{i_v}+gp_inc)/(rho_f*omega^2);
                    else
                        displacement{i_v} = gScalarField_p{i_v}/(rho_f*omega^2);
                    end
                end
            else
                totField = scalarField{i_v};
                if d_p == 3
                    rho_f = varCol{i_v}.rho; 
                    displacement{i_v} = gScalarField_p{i_v}/(rho_f*omega^2);  
                end
            end
            if analyticSolutionExist
                p = layer{i_v}.p;
                data.analytic = real(makeDynamic(p, para{i_v}, omega)); 
                p_e = p-scalarField{i_v};
                p2 = abs(p).^2;
                p_e2 = abs(p_e).^2;
                data.Error = sqrt(p_e2/max(p2));
                if d_p == 3
                    dp = [layer{i_v}.dpdx, layer{i_v}.dpdy, layer{i_v}.dpdz];
                    dp_e = dp-gScalarField_p{i_v};

                    p2 = abs(p).^2;
                    dp2 = sum(abs(dp).^2,2);
                    p_e2 = abs(p_e).^2;
                    dp_e2 = sum(abs(dp_e).^2,2);

                    k = varCol{i_v}.k;
                    data.Error = sqrt(p_e2/max(p2));
                    data.ErrorGrad = sqrt(dp_e2/max(dp2));
                    data.ErrorEnergy = sqrt((dp_e2 + k^2*p_e2)/max(dp2 + k^2*p2));
                end
            end
            data.totField = real(makeDynamic(totField, para{i_v}, omega));
            data.totFieldAbs = abs(makeDynamic(totField, para{i_v}, omega));

        case 'solid'
            if analyticSolutionExist
                u = [layer{2}.u_x,layer{2}.u_y,layer{2}.u_z];                
                u_e = u-displacement{i_v};

                u2 = sum(abs(u).^2,2);
                u_e2 = sum(abs(u_e).^2,2);

                sigma = [layer{2}.sigma_xx,layer{2}.sigma_yy,layer{2}.sigma_zz,layer{2}.sigma_yz,layer{2}.sigma_xz,layer{2}.sigma_xy];
                C = varCol{i_v}.C;

                strain_vec = (C\sigma.').';

                strain_e = strain_vec-strain_h{i_v};

                eCe = real(sum((strain_e*C).*conj(strain_e),2)); % the usage of real() is to remove machine epsilon imaginary part
                uCu = real(sum((strain_vec*C).*conj(strain_vec),2)); % the usage of real() is to remove machine epsilon imaginary part

                rho_s = varCol{i_v}.rho; 
                data.Error = sqrt(u_e2/max(u2));
                data.ErrorGrad = sqrt(eCe/max(uCu));
                data.ErrorEnergy = sqrt((eCe + rho_s*omega^2*u_e2)/max(uCu + rho_s*omega^2*u2));

                data.analytic = real(makeDynamic(u, para{i_v}, omega)); 
            end
            data.stress = real(makeDynamic(strain_h{i_v}*C, para{i_v}, omega));
            if withDensity
                data.density = real(makeDynamic(density{i_v}, para{i_v}, omega));
            end
    end
    fprintf('using %12f seconds.', toc)
    data.displacement = real(makeDynamic(displacement{i_v}, para{i_v}, omega));
    
    fprintf(['\n%-' num2str(stringShift) 's'], '    Creating VTK-file ... ')
    para{i_v}.name = [para{i_v}.name '_' num2str(i_v)];
    tic
    makeVTKfile(data, para{i_v});
    fprintf('using %12f seconds.', toc)
end

if options.para_options.plotMesh
    fprintf(['\n%-' num2str(stringShift) 's'], '    Creating mesh-files ... ')
    tic
    createVTKmeshFiles(varCol, 'para_options', options.para_options)
    fprintf('using %12f seconds.', toc)
end
