function [nodes, noNodes, visElements, cornerNode, ...
          noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, solid)

XiVec = insertUniform(Xi, extraXiPts);
EtaVec = insertUniform(Eta, extraEtaPts);
ZetaVec = insertUniform(Zeta, extraZetaPts);

% number of distinct knot values

noXiKnots = length(XiVec);
noEtaKnots = length(EtaVec);
noZetaKnots = length(ZetaVec);

noNodes = noXiKnots*noEtaKnots*noZetaKnots;
nodes  = zeros(noNodes,3);
% nodes  = zeros(noNodes,4);
cornerNode = zeros(noElems,1);

count = 1;
e = 1;
for k = 1:noZetaKnots
    zeta = ZetaVec(k);
    for j = 1:noEtaKnots
        eta = EtaVec(j);
        for i = 1:noXiKnots
            xi = XiVec(i);

            v = evaluateNURBS(solid, [xi, eta, zeta]);

            nodes(count,1) = v(1);
            nodes(count,2) = v(2);
            nodes(count,3) = v(3);
%             nodes(count,4) = count;
            
            if ismember(xi, Xi) && ismember(eta, Eta) && ismember(zeta, Zeta) && ...
                  xi ~= Xi(end) && eta ~= Eta(end) && zeta ~= Zeta(end)
                cornerNode(e) = count;
                e = e+1;
            end
            count = count + 1;
        end
    end
end

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



