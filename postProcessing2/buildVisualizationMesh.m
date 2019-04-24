function visObj = buildVisualizationMesh(nurbs, extraPts)

switch nurbs.type
    case '3Dsurface'
        
        XiVec = insertUniform(nurbs.knots{1}, extraPts(1));
        EtaVec = insertUniform(nurbs.knots{2}, extraPts(2));

        noXiKnots = length(XiVec);
        noEtaKnots = length(EtaVec);

        noNodes = noXiKnots*noEtaKnots;
        nodes_x = zeros(noXiKnots,noEtaKnots);
        nodes_y = zeros(noXiKnots,noEtaKnots);
        nodes_z = zeros(noXiKnots,noEtaKnots);

        parfor j = 1:noEtaKnots
            nodes_x_temp = zeros(noXiKnots,1);
            nodes_y_temp = zeros(noXiKnots,1);
            nodes_z_temp = zeros(noXiKnots,1);
            eta = EtaVec(j);
            for i = 1:noXiKnots
                xi = XiVec(i);

                v = evaluateNURBS(nurbs, [xi, eta, 1]);

                nodes_x_temp(i) = v(1);
                nodes_y_temp(i) = v(2);
                nodes_z_temp(i) = v(3);
            end
            nodes_x(:,j) = nodes_x_temp;
            nodes_y(:,j) = nodes_y_temp;
            nodes_z(:,j) = nodes_z_temp;
        end
        nodes = [reshape(nodes_x, noNodes, 1), reshape(nodes_y, noNodes, 1), reshape(nodes_z, noNodes, 1)];


        noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
        visElements = zeros(noVisElems,4);
        eVis = 1;

        for j = 1:noEtaKnots-1
            for i = 1:noXiKnots-1
                visElements(eVis,1) = i   +   (j-1)*noXiKnots;
                visElements(eVis,2) = i+1 +   (j-1)*noXiKnots;
                visElements(eVis,3) = i+1 +       j*noXiKnots;
                visElements(eVis,4) = i   +       j*noXiKnots;

                eVis = eVis + 1;
            end
        end

        
    case '3Dvolume'
        uniqueZeta = unique(nurbs.knots{3});
        
        XiVec = insertUniform(nurbs.knots{1}, extraPts(1));
        EtaVec = insertUniform(nurbs.knots{2}, extraPts(2));
        ZetaVec = insertUniform(nurbs.knots{3}, extraPts(3));

        noXiKnots = length(XiVec);
        noEtaKnots = length(EtaVec);
        noZetaKnots = length(ZetaVec);

        noNodes = noXiKnots*noEtaKnots*noZetaKnots;
        nodes_x  = zeros(noXiKnots*noEtaKnots,noZetaKnots);
        nodes_y  = zeros(noXiKnots*noEtaKnots,noZetaKnots);
        nodes_z  = zeros(noXiKnots*noEtaKnots,noZetaKnots);

        parfor k = 1:noZetaKnots
            zeta = ZetaVec(k);
            nodes_x_temp = zeros(noEtaKnots*noXiKnots,1);
            nodes_y_temp = zeros(noEtaKnots*noXiKnots,1);
            nodes_z_temp = zeros(noEtaKnots*noXiKnots,1);
            count = 1;
            for j = 1:noEtaKnots
                eta = EtaVec(j);
                for i = 1:noXiKnots
                    xi = XiVec(i);

                    v = evaluateNURBS(nurbs, [xi, eta, zeta]);

                    nodes_x_temp(count) = v(1);
                    nodes_y_temp(count) = v(2);
                    nodes_z_temp(count) = v(3);

                    count = count + 1;
                end
            end
            nodes_x(:,k) = nodes_x_temp;
            nodes_y(:,k) = nodes_y_temp;
            nodes_z(:,k) = nodes_z_temp;
        end
        nodes = [reshape(nodes_x, noNodes, 1), reshape(nodes_y, noNodes, 1), reshape(nodes_z, noNodes, 1)];

        noVisElems  = (noXiKnots-1)*(noEtaKnots-1)*(noZetaKnots-1);
        visElements = zeros(noVisElems,8);
        eVis = 1;
        for k = 1:noZetaKnots-1
            for j = 1:noEtaKnots-1
                for i = 1:noXiKnots-1
                    visElements(eVis,1) = i   +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                    visElements(eVis,2) = i+1 +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                    visElements(eVis,3) = i+1 +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                    visElements(eVis,4) = i   +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
                    visElements(eVis,5) = i   +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
                    visElements(eVis,6) = i+1 +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
                    visElements(eVis,7) = i+1 +       j*noXiKnots +       k*noEtaKnots*noXiKnots;
                    visElements(eVis,8) = i   +       j*noXiKnots +       k*noEtaKnots*noXiKnots;

                    eVis = eVis + 1;
                end
            end
        end
        visObj.ZetaVec = ZetaVec;
        visObj.noZetaKnots = noZetaKnots;
        
        nurbs2 = nurbs;
        nurbs2.knots{1} = [XiVec(1), XiVec.', XiVec(end)];
        nurbs2.knots{2} = [EtaVec(1), EtaVec.', EtaVec(end)];
        nurbs2.knots{3} = [ZetaVec(1), ZetaVec.', ZetaVec(end)];
        nurbs2.number = [length(XiVec),length(EtaVec),length(ZetaVec)];
        nurbs2.degree = [1,1,1];
        varCol.nurbs = nurbs2;
        varCol = generateIGA3DMesh_new(varCol);
        
        visObj.index = varCol.index;
        visObj.noElems = varCol.noElems;
        visObj.elRangeXi = varCol.elRangeXi;
        visObj.elRangeEta = varCol.elRangeEta;
        visObj.elRangeZeta = varCol.elRangeZeta;
end
visObj.nodes = nodes;
visObj.noNodes = noNodes;
visObj.visElements = visElements;
visObj.noXiKnots = noXiKnots;
visObj.noEtaKnots = noEtaKnots;
visObj.XiVec = XiVec;
visObj.EtaVec = EtaVec;


