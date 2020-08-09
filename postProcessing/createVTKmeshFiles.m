function createVTKmeshFiles(varargin)
varCol = varargin{1};
options = struct('U', NaN, ... 
                 'para_options', []);
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
    U = options.U{i_v};

    extraXiPts = para.extraXiPts;
    extraEtaPts = para.extraEtaPts;
    extraZetaPts = para.extraZetaPts;
    degree = varCol{i_v}.degree; % assume p_xi is equal in all patches
    Eps = 1e4*eps;
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
    if isfield('varCol{i_v}','omega')
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
    if strcmp(varCol{i_v}.media,'fluid') && para.plotDisplacementVectors
        warning('displacement (based on grad(p)) is not yet implemented for fluids meshes')
        plotDisplacementVectors = false;
    else
        plotDisplacementVectors = para.plotDisplacementVectors;
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

            options = struct('name',[para.name '_' num2str(i_v) '_mesh'], 'celltype', 'VTK_POLY_LINE');

            data.nodes = nodes;
            data.visElements = visElements;

            makeVTKfile(data, options);
        case 3
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
                        container_e{counter}.displacement = R{1}*Usctr;
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
                        container_e{counter}.displacement = R{1}*Usctr;
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
                        container_e{counter}.displacement = R{1}*Usctr;
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
                    if d == 3
                        displacement(nodesCount+1:nodesCount+container{e}.container_e{i}.noNodes,:) = container{e}.container_e{i}.displacement;
                    end
                    nodesCount = nodesCount + container{e}.container_e{i}.noNodes;
                end
            end

            para = struct('name',[para.name '_' num2str(i_v) '_mesh'], 'celltype', 'VTK_POLY_LINE', ...
                'plotDisplacementVectors', plotDisplacementVectors,'plotTimeOscillation', para.plotTimeOscillation);

            data.nodes = nodes;
            data.displacement = real(makeDynamic(displacement, para, omega));
            data.visElements = visElements;

            makeVTKfile(data, para);
    end
end