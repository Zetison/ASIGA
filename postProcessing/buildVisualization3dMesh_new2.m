function [nodes, noNodes, visElements, cornerNode, ...
          noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new2(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, solid)

XiVec = insertUniform(Xi, extraXiPts);
EtaVec = insertUniform(Eta, extraEtaPts);
ZetaVec = insertUniform(Zeta, extraZetaPts);

% number of distinct knot values

noXiKnots = length(XiVec);
noEtaKnots = length(EtaVec);
noZetaKnots = length(ZetaVec);

noNodes = noXiKnots*noEtaKnots*noZetaKnots;
nodes  = zeros(noNodes,4);
cornerNode = zeros(noElems,1);

count = 1;
e = 1;
chan = zeros(noXiKnots, noEtaKnots, noZetaKnots);
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
            nodes(count,4) = count;
            
            if ismember(xi, Xi) && ismember(eta, Eta) && ismember(zeta, Zeta) && ...
                  xi ~= Xi(end) && eta ~= Eta(end) && zeta ~= Zeta(end)
                cornerNode(e) = count;
                e = e+1;
            end
            chan(i,j,k) = count;
            count = count + 1;
        end
    end
end
noElemsXi = length(unique(Xi))-1;
noElemsEta = length(unique(Eta))-1;
noElemsZeta = length(unique(Zeta))-1;

elConnXi = buildConnectivityVisualization(extraXiPts, noElemsXi);
elConnEta = buildConnectivityVisualization(extraEtaPts, noElemsEta);
elConnZeta = buildConnectivityVisualization(extraZetaPts, noElemsZeta);

e = 1;
noVisElemsPerElems = (extraXiPts+1)*(extraEtaPts+1)*(extraZetaPts+1);
visElements = zeros(noElems, noVisElemsPerElems,8);
for kk = 1:noElemsZeta
    zetaConn = elConnZeta(kk,:);
    for jj = 1:noElemsEta
        etaConn = elConnEta(jj,:);
        for ii = 1:noElemsXi
            xiConn = elConnXi(ii,:);
            
            eVis = 1;            
            for k = 1:extraZetaPts+1
                for j = 1:extraEtaPts+1
                    for i = 1:extraXiPts+1
                        cornerIdx = chan(xiConn(i), etaConn(j), zetaConn(k));
                        visElements(e,eVis,1) = cornerIdx;
                        visElements(e,eVis,2) = cornerIdx + 1;
                        visElements(e,eVis,3) = cornerIdx + 1  + noXiKnots;
                        visElements(e,eVis,4) = cornerIdx      + noXiKnots;
                        visElements(e,eVis,5) = cornerIdx                    + noEtaKnots*noXiKnots;
                        visElements(e,eVis,6) = cornerIdx + 1                + noEtaKnots*noXiKnots;
                        visElements(e,eVis,7) = cornerIdx + 1  + noXiKnots   + noEtaKnots*noXiKnots;
                        visElements(e,eVis,8) = cornerIdx      + noXiKnots   + noEtaKnots*noXiKnots;
                        
                        eVis = eVis + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end
% 
% 
% % noVisElems  = (noXiKnots-1)*(noEtaKnots-1)*(noZetaKnots-1);
% keyboard
% for e = 1:noElems
%     for eVis = 1:size(elementToGlobInd,2)
%         indices = elementToGlobInd(e,eVis,:);
%         i = indices(1);
%         j = indices(2);
%         k = indices(3);
%         visElements(e,eVis,1) = i   +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
%         visElements(e,eVis,2) = i+1 +   (j-1)*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
%         visElements(e,eVis,3) = i+1 +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
%         visElements(e,eVis,4) = i   +       j*noXiKnots +   (k-1)*noEtaKnots*noXiKnots;
%         visElements(e,eVis,5) = i   +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
%         visElements(e,eVis,6) = i+1 +   (j-1)*noXiKnots +       k*noEtaKnots*noXiKnots;
%         visElements(e,eVis,7) = i+1 +       j*noXiKnots +       k*noEtaKnots*noXiKnots;
%         visElements(e,eVis,8) = i   +       j*noXiKnots +       k*noEtaKnots*noXiKnots;
%     end
% end
