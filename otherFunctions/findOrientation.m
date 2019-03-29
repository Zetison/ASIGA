function rightHandedOrientation = findOrientation(varCol)

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);

index = varCol.index;

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;

e = 1;
pts = controlPts(element(e,:),:);
idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
idEta = index(e,2);
idZeta = index(e,3);

% calculate element ranges
% elRangeU, elRangeV and elRangeW comes from generateIGA3DMesh
Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]
xi   = parent2ParametricSpace(Xi_e,  0.5);
eta  = parent2ParametricSpace(Eta_e, 0.5);
zeta = parent2ParametricSpace(Zeta_e,0.5);
[~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
J_1 = det(pts'*[dRdxi' dRdeta' dRdzeta']);
if J_1 > 0
    rightHandedOrientation = true;
else
    rightHandedOrientation = false;
end