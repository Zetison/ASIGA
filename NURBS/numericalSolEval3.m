function [u, v] = numericalSolEval3(xi, eta, p, q, Xi, Eta, controlPts, weights, Ux, Uy, Uz, elRangeXi, elRangeEta, noElemsXi, XiEtaMesh,zeta1Nodes)


[R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

idx_xi = findElementInMesh(xi, elRangeXi);
idx_eta = findElementInMesh(eta, elRangeEta);

e = noElemsXi*(idx_eta-1) + idx_xi;

sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector
v = R_fun*controlPts(sctrXiEta,:);
u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];