function createVTKmeshFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options)

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
type = patches{1}.nurbs.type;
d = varCol.dimension;
omega = varCol.omega;
isOuterDomain = varCol.isOuterDomain;
if isOuterDomain
    model = varCol.model;
end
switch type
    case '3Dsurface'
        %% Create Mesh vtk files
        container = cell(1,noElems);
%         for e = 1:noElems
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
            noVisElems  = 2*noXiKnots+2*(noEtaKnots-2)+1;
            visElements_e = 1:noVisElems;
            
            noNodes = noVisElems;    
            xiKnots = linspace(Xi_e(1),Xi_e(2)-eps,noXiKnots);        
            etaKnots = linspace(Eta_e(1),Eta_e(2)-eps,noEtaKnots);
            counter = 1;
            Eps = 1e10*eps;
            Eps2 = 1e13*eps;
            nodes_e = zeros(noNodes,3);
            for i = 1:noXiKnots
                xi = xiKnots(i);
                [~, v] = numericalSolEval_final_surf(xi, etaKnots(1), p_xi, p_eta, Xi, Eta, wgts, pts, Usctr);
                [~, dRdxi, dRdeta] = NURBS2DBasis(min(max(xi-Eps,xiKnots(1)+Eps),xiKnots(end)-Eps), etaKnots(1)+10*eps, p_xi, p_eta, Xi, Eta, wgts);
                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd/J_1;
                if any(isnan(n)) || any(isinf(n))
                    error('The normal vector contains nan or inf values!')
                end
                nodes_e(counter,:) = v+n*Eps2;
                counter = counter + 1;
            end
            for j = 2:noEtaKnots
                eta = etaKnots(j);
                [~, v] = numericalSolEval_final_surf(xiKnots(end), eta, p_xi, p_eta, Xi, Eta, wgts, pts, Usctr);
                [~, dRdxi, dRdeta] = NURBS2DBasis(xiKnots(end)-Eps, eta-Eps, p_xi, p_eta, Xi, Eta, wgts);
                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd/J_1;
                if any(isnan(n)) || any(isinf(n))
                    error('The normal vector contains nan or inf values!')
                end
                nodes_e(counter,:) = v+n*Eps2;
                counter = counter + 1;
            end
            for i = (noXiKnots-1):-1:1
                xi = xiKnots(i);
                [~, v] = numericalSolEval_final_surf(xi, etaKnots(end), p_xi, p_eta, Xi, Eta, wgts, pts, Usctr);
                [~, dRdxi, dRdeta] = NURBS2DBasis(xi+Eps, etaKnots(end)-Eps, p_xi, p_eta, Xi, Eta, wgts);
                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd/J_1;
                if any(isnan(n)) || any(isinf(n))
                    error('The normal vector contains nan or inf values!')
                end
                nodes_e(counter,:) = v+n*Eps2;
                counter = counter + 1;
            end
            for j = (noEtaKnots-1):-1:1
                eta = etaKnots(j);
                [~, v] = numericalSolEval_final_surf(xiKnots(1), eta, p_xi, p_eta, Xi, Eta, wgts, pts, Usctr);
                [~, dRdxi, dRdeta] = NURBS2DBasis(xiKnots(1)+Eps, eta+Eps, p_xi, p_eta, Xi, Eta, wgts);
                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);
                n = crossProd/J_1;
                if any(isnan(n)) || any(isinf(n))
                    error('The normal vector contains nan or inf values!')
                end
                nodes_e(counter,:) = v+n*Eps2;
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

        options = struct('name',[options.name 'mesh_'], 'celltype', 'VTK_POLY_LINE');
            
        data.nodes = nodes;
        data.visElements = visElements;


        makeVTKfile(data, options);
    case '3Dvolume'
        Zeta = nurbs.knots{3};
        %% Create Mesh vtk files
        noXiKnots = visObj.noXiKnots;
        noEtaKnots = visObj.noEtaKnots;
        noZetaKnots = visObj.noZetaKnots;
        uniqueXiKnots = unique(Xi);
        uniqueEtaKnots = unique(Eta);
        uniqueZetaKnots = unique(Zeta);
        noUniqueXiKnots = numel(uniqueXiKnots);
        noUniqueEtaKnots = numel(uniqueEtaKnots);
        noUniqueZetaKnots = numel(uniqueZetaKnots);
        meshNodes = zeros(noUniqueZetaKnots*noUniqueEtaKnots*noXiKnots ...
                         +noUniqueZetaKnots*noEtaKnots*noUniqueXiKnots ...
                         +noZetaKnots*noUniqueEtaKnots*noUniqueXiKnots,3);
        if exist('U','var')
            displacement = zeros(noUniqueZetaKnots*noUniqueEtaKnots*noXiKnots ...
                             +noUniqueZetaKnots*noEtaKnots*noUniqueXiKnots ...
                             +noZetaKnots*noUniqueEtaKnots*noUniqueXiKnots,3);
        end
        count = 1;
        count2 = 1;
        visElements = cell(noUniqueEtaKnots*noUniqueZetaKnots + noUniqueXiKnots*noUniqueZetaKnots + noUniqueXiKnots*noUniqueEtaKnots, 1);
        % const eta and zeta
        for k = 1:noUniqueZetaKnots
            for j = 1:noUniqueEtaKnots
                visElements{count2} = count:(count+noXiKnots-1);
                count2 = count2 + 1;
                for i = 1:noXiKnots
                    zeta = uniqueZetaKnots(k);
                    eta = uniqueEtaKnots(j);
                    xi = XiVec(i);    

                    if exist('U','var')
                        [u, v] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                        if d == 3
                            displacement(count,:) = u.';
                        end
                    else
                        v = evaluateNURBS(nurbs,[xi,eta,zeta]);
                    end


                    meshNodes(count,:) = v';

                    count = count + 1;
                end
            end
        end

        % const xi and zeta
        for k = 1:noUniqueZetaKnots
            for i = 1:noUniqueXiKnots 
                visElements{count2} = count:(count+noEtaKnots-1);
                count2 = count2 + 1;
                for j = 1:noEtaKnots
                    zeta = uniqueZetaKnots(k);
                    eta = EtaVec(j);
                    xi = uniqueXiKnots(i);    

                    if exist('U','var')
                        [u, v] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                        if d == 3
                            displacement(count,:) = u.';
                        end
                    else
                        v = evaluateNURBS(nurbs,[xi,eta,zeta]);
                    end

                    meshNodes(count,:) = v';

                    count = count + 1;
                end
            end
        end

        % const xi and eta 
        for j = 1:noUniqueEtaKnots
            for i = 1:noUniqueXiKnots   
                visElements{count2} = count:(count+noZetaKnots-1);
                count2 = count2 + 1;           
                for k = 1:noZetaKnots
                    zeta = ZetaVec(k);
                    eta = uniqueEtaKnots(j);
                    xi = uniqueXiKnots(i);       

                    if exist('U','var')
                        [u, v] = numericalSolEval_final(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights, controlPts, U);
                        if d == 3
                            displacement(count,:) = u.';
                        end
                    else
                        v = evaluateNURBS(nurbs,[xi,eta,zeta]);
                    end

                    meshNodes(count,:) = v';

                    count = count + 1;
                end
            end
        end

        if exist('U','var')
            if d == 1
                options = struct('name',[options.name 'mesh_'], 'celltype', 'VTK_POLY_LINE',  'plotTimeOscillation', options.plotTimeOscillation);
            else
                options = struct('name',[options.name 'mesh_'], 'celltype', 'VTK_POLY_LINE', 'plotTimeOscillation', options.plotTimeOscillation, 'plotDisplacementVectors',1);
            end
            if max(max(abs(displacement))) < 1e-6
                data.displacement = 0.25*displacement*4e9;
            else
                data.displacement = displacement;
            end
        else
            options = struct('name',[options.name 'mesh_'], 'celltype', 'VTK_POLY_LINE');
        end

        data.nodes = meshNodes;
        data.visElements = visElements;


        makeVTKfile(data, options);
end