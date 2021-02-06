function createVTKmeshFiles(varargin)
varCol = varargin{1};
options = struct('para_options', []);
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
for i_v = 1:numel(varCol)
    para = options.para_options;
    U = varCol{i_v}.U;
    isOuterDomain = i_v == 1;

    extraXiPts = para.extraXiPts;
    extraEtaPts = para.extraEtaPts;
    extraZetaPts = para.extraZetaPts;
    degree = varCol{i_v}.degree; % assume p_xi is equal in all patches
    n_en = prod(degree+1);
    
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
        Ux = U(1:d:noDofs);
        Uy = U(2:d:noDofs);
        Uz = U(3:d:noDofs);
        U = [Ux, Uy, Uz];
    end
    d_p = patches{1}.nurbs.d_p;
    switch d_p
        case 2
            container = cell(1,noElems);
    %         for e = 1:noElems
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
                
                noXiKnots = 2+extraXiPts;
                noEtaKnots = 2+extraEtaPts;
                noVisElems  = 2*noXiKnots+2*(noEtaKnots-2)+1;
                visElements_e = 1:noVisElems;

                noNodes = noVisElems;    
                xiKnots = linspace(Xi_e(1,1),Xi_e(1,2)-eps,noXiKnots);        
                etaKnots = linspace(Xi_e(2,1),Xi_e(2,2)-eps,noEtaKnots);
                counter = 1;
                nodes_e = zeros(noNodes,3);
                I = findKnotSpans(degree, [xiKnots(1),etaKnots(1)], knots);

                for i = 1:noXiKnots
                    xi = xiKnots(i);
                    R = NURBSbasis(I,[min(max(xi,xiKnots(1)),xiKnots(end)), etaKnots(1)], degree, knots, wgts);
                    nodes_e(counter,:) = R{1}*pts;
                    counter = counter + 1;
                end
                for j = 2:noEtaKnots
                    eta = etaKnots(j);
                    R = NURBSbasis(I,[xiKnots(end), eta], degree, knots, wgts);
                    nodes_e(counter,:) = R{1}*pts;
                    counter = counter + 1;
                end
                for i = (noXiKnots-1):-1:1
                    xi = xiKnots(i);
                    R = NURBSbasis(I,[xi, etaKnots(end)], degree, knots, wgts);
                    nodes_e(counter,:) = R{1}*pts;
                    counter = counter + 1;
                end
                for j = (noEtaKnots-1):-1:1
                    eta = etaKnots(j);
                    R = NURBSbasis(I,[xiKnots(1), eta], degree, knots, wgts);
                    nodes_e(counter,:) = R{1}*pts;
                    counter = counter + 1;
                end
                container{e}.nodes = nodes_e;
                container{e}.noNodes = noNodes;
                container{e}.noVisElems = noVisElems;
                container{e}.visElements = visElements_e;
            end
            noNodes = 0;
            for e = 1:noElems
                noNodes = noNodes + container{e}.noNodes;
            end
            visElements = cell(noElems,1);
            nodes = zeros(noNodes,3);
            nodesCount = 0;
            for e = 1:noElems
                visElements{e} = nodesCount + container{e}.visElements;
                nodes(nodesCount+1:nodesCount+container{e}.noNodes,:) = container{e}.nodes;
                nodesCount = nodesCount + container{e}.noNodes;
            end

            data.nodes = nodes;
            data.visElements = visElements;

            makeVTKfile(data, 'name', [para.name '_' num2str(i_v) '_mesh'], 'celltype', 'VTK_POLY_LINE');
        case 3
            Eps = 1e-6;
            container = cell(1,noElems);
%             for e = 1:noElems
            parfor e = 1:noElems
                patch = pIndex(e);
                knots = knotVecs{patch};

                Xi_e = zeros(d_p,2);
                for i = 1:d_p
                    Xi_e(i,:) = elRange{i}(index(e,i),:);
                end
                Xi_e(1,2) = Xi_e(1,2)-eps;
                Xi_e(2,2) = Xi_e(2,2)-eps;
                Xi_e(3,2) = Xi_e(3,2)-eps;

                sctr = element(e,:);
                pts = controlPts(sctr,:);
                wgts = weights(element2(e,:),:); % New
                Usctr = U(sctr,:);

                noXiKnots = 2+extraXiPts;
                noEtaKnots = 2+extraEtaPts;
                noZetaKnots = 2+extraZetaPts;

                xi = linspace(Xi_e(1,1),Xi_e(1,2)-eps,noXiKnots).';        
                eta = linspace(Xi_e(2,1),Xi_e(2,2)-eps,noEtaKnots).';
                zeta = linspace(Xi_e(3,1),Xi_e(3,2)-eps,noZetaKnots).';
                I = findKnotSpans(degree, [xi(1),eta(1),zeta(1)], knots);
                container_e = cell(12,1);
                counter = 1;
                for i = 1:2
                    for j = 1:2
                        R = NURBSbasis(I,[Xi_e(1,i)*ones(noZetaKnots,1), Xi_e(2,j)*ones(noZetaKnots,1), zeta], degree, knots, wgts);
                        container_e{counter}.nodes = R{1}*pts;
                        if d == 3
                            container_e{counter}.displacement = R{1}*Usctr;
                        else
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
                            
                            dudx = dRdx*Usctr;
                            dudy = dRdy*Usctr;
                            dudz = dRdz*Usctr;
                            singularJ = abs(J_1) < Eps;
                            dudx(or(singularJ, isnan(dudx))) = 0;
                            dudy(or(singularJ, isnan(dudy))) = 0;
                            dudz(or(singularJ, isnan(dudz))) = 0;
                            container_e{counter}.gScalarField_p = [dudx,dudy,dudz];
                        end
                        container_e{counter}.noNodes = noZetaKnots;
                        container_e{counter}.visElements = 1:noZetaKnots;
                        container_e{counter}.noVisElems = noZetaKnots;
                        counter = counter + 1;
                    end
                end
                for i = 1:2
                    for j = 1:2
                        R = NURBSbasis(I,[Xi_e(1,i)*ones(noEtaKnots,1), eta, Xi_e(3,j)*ones(noEtaKnots,1)], degree, knots, wgts);
                        container_e{counter}.nodes = R{1}*pts;
                        if d == 3
                            container_e{counter}.displacement = R{1}*Usctr;
                        else
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
                            
                            dudx = dRdx*Usctr;
                            dudy = dRdy*Usctr;
                            dudz = dRdz*Usctr;
                            singularJ = abs(J_1) < Eps;
                            dudx(or(singularJ, isnan(dudx))) = 0;
                            dudy(or(singularJ, isnan(dudy))) = 0;
                            dudz(or(singularJ, isnan(dudz))) = 0;
                            container_e{counter}.gScalarField_p = [dudx,dudy,dudz];
                        end
                        container_e{counter}.noNodes = noEtaKnots;
                        container_e{counter}.visElements = 1:noEtaKnots;
                        container_e{counter}.noVisElems = noEtaKnots;
                        counter = counter + 1;
                    end
                end
                for i = 1:2
                    for j = 1:2
                        R = NURBSbasis(I,[xi, Xi_e(2,i)*ones(noXiKnots,1), Xi_e(3,j)*ones(noXiKnots,1)], degree, knots, wgts);
                        
                        container_e{counter}.nodes = R{1}*pts;
                        if d == 3
                            container_e{counter}.displacement = R{1}*Usctr;
                        else
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
                            
                            dudx = dRdx*Usctr;
                            dudy = dRdy*Usctr;
                            dudz = dRdz*Usctr;
                            singularJ = abs(J_1) < Eps;
                            dudx(or(singularJ, isnan(dudx))) = 0;
                            dudy(or(singularJ, isnan(dudy))) = 0;
                            dudz(or(singularJ, isnan(dudz))) = 0;
                            container_e{counter}.gScalarField_p = [dudx,dudy,dudz];
                        end
                        container_e{counter}.noNodes = noXiKnots;
                        container_e{counter}.visElements = 1:noXiKnots;
                        container_e{counter}.noVisElems = noXiKnots;
                        counter = counter + 1;
                    end
                end
                container{e}.container_e = container_e;
            end
            noNodes = 0;
            for e = 1:noElems
                for i = 1:12
                    noNodes = noNodes + container{e}.container_e{i}.noNodes;
                end
            end
            visElements = cell(noElems*12,1);
            nodes = zeros(noNodes,3);
            displacement = zeros(noNodes,3);
            nodesCount = 0;
            for e = 1:noElems
                for i = 1:12
                    visElements{12*(e-1)+i} = nodesCount + container{e}.container_e{i}.visElements;
                    nodes(nodesCount+1:nodesCount+container{e}.container_e{i}.noNodes,:) = container{e}.container_e{i}.nodes;
                    switch varCol{i_v}.media
                        case 'fluid'
                            displacement(nodesCount+1:nodesCount+container{e}.container_e{i}.noNodes,:) = container{e}.container_e{i}.gScalarField_p;
                        case 'solid'
                            displacement(nodesCount+1:nodesCount+container{e}.container_e{i}.noNodes,:) = container{e}.container_e{i}.displacement;
                    end
                    nodesCount = nodesCount + container{e}.container_e{i}.noNodes;
                end
            end
            switch varCol{i_v}.media
                case 'fluid'
                    if isOuterDomain
                        if d_p == 3
                            rho_f = varCol{i_v}.rho;
                            if varCol{1}.splitExteriorFields
                                gp_inc = [varCol{i_v}.dp_incdx_(nodes),varCol{i_v}.dp_incdy_(nodes),varCol{i_v}.dp_incdz_(nodes)];
                                displacement = (displacement+gp_inc)/(rho_f*omega^2);
                            else
                                displacement = displacement/(rho_f*omega^2);
                            end
                        end
                    else
                        if d_p == 3
                            rho_f = varCol{i_v}.rho; 
                            displacement = displacement/(rho_f*omega^2);  
                        end
                    end
            end
            
            data.nodes = nodes;
            data.displacement = real(makeDynamic(displacement, para, omega));
            data.visElements = visElements;

            makeVTKfile(data, 'name',[para.name '_' num2str(i_v) '_mesh'], 'celltype', 'VTK_POLY_LINE', ...
                'plotDisplacementVectors', para.plotDisplacementVectors,'plotTimeOscillation', para.plotTimeOscillation);
    end
end