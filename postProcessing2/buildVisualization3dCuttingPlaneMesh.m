function [nodes, noNodes, visElements, ...
          noEtaKnots, noZetaKnots, ...
          EtaVec, ZetaVec] = buildVisualization3dCuttingPlaneMesh(Eta, Zeta, extraEtaPts, extraZetaPts, nurbs)

EtaVec = insertUniform(Eta, extraEtaPts);
ZetaVec = insertUniform(Zeta, extraZetaPts);

% number of distinct knot values

noEtaKnots = length(EtaVec);
noZetaKnots = length(ZetaVec);

noNodes = noEtaKnots*noZetaKnots;
nodes_x = zeros(noEtaKnots,noZetaKnots);
nodes_y = zeros(noEtaKnots,noZetaKnots);
nodes_z = zeros(noEtaKnots,noZetaKnots);
nodes = [];
for xi = [0, 0.5]
    parfor k = 1:noZetaKnots
        nodes_x_temp = zeros(noEtaKnots,1);
        nodes_y_temp = zeros(noEtaKnots,1);
        nodes_z_temp = zeros(noEtaKnots,1);
        zeta = ZetaVec(k);
        for j = 1:noEtaKnots
            eta = EtaVec(j);

            v = evaluateNURBS(nurbs, [xi, eta, zeta]);

            nodes_x_temp(j) = v(1);
            nodes_y_temp(j) = v(2);
            nodes_z_temp(j) = v(3);
        end
        nodes_x(:,k) = nodes_x_temp;
        nodes_y(:,k) = nodes_y_temp;
        nodes_z(:,k) = nodes_z_temp;
    end
    nodes = [nodes; reshape(nodes_x, noNodes, 1), reshape(nodes_y, noNodes, 1), reshape(nodes_z, noNodes, 1)];
end


noVisElems  = (noEtaKnots-1)*(noZetaKnots-1);
visElements = zeros(2*noVisElems,4);
eVis = 1;

for i = [0 1]
    for k = 1:noZetaKnots-1
        for j = 1:noEtaKnots-1
            visElements(eVis,1) = j   +   (k-1)*noEtaKnots + i*noNodes;
            visElements(eVis,2) = j+1 +   (k-1)*noEtaKnots + i*noNodes;
            visElements(eVis,3) = j+1 +       k*noEtaKnots + i*noNodes;
            visElements(eVis,4) = j   +       k*noEtaKnots + i*noNodes;

            eVis = eVis + 1;
        end
    end
end


