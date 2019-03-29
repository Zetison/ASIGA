function [nodes, noNodes, visElements, ...
          noXiKnots, noEtaKnots, ...
          XiVec, EtaVec] = buildVisualization3dSurfaceMesh(Xi, Eta, extraXiPts, extraEtaPts, nurbs)

XiVec = insertUniform(Xi, extraXiPts);
EtaVec = insertUniform(Eta, extraEtaPts);

% number of distinct knot values

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



