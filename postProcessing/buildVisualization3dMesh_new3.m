function [nodes, noNodes, visElements, ...
          noXiKnots, noEtaKnots, noZetaKnots, ...
          XiVec, EtaVec, ZetaVec] = buildVisualization3dMesh_new3(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, nurbs)

XiVec = insertUniform(Xi, extraXiPts);
EtaVec = insertUniform(Eta, extraEtaPts);
ZetaVec = insertUniform(Zeta, extraZetaPts);

% number of distinct knot values

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



