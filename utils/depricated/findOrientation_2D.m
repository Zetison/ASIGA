function rightHandedOrientation = findOrientation_2D(varCol)
error('Depricated, ASIGA assumes right handed orientation')
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

index = varCol.index;

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;

e = 1;
pts = controlPts(element(e,:),:);
idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
idEta = index(e,2);

% calculate element ranges
% elRangeU, elRangeV and elRangeW comes from generateIGA3DMesh
Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]
xi   = parent2ParametricSpace(Xi_e,  0.5);
eta  = parent2ParametricSpace(Eta_e, 0.5);
[~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
J_1 = det(pts'*[dRdxi' dRdeta']);
if J_1 > 0
    rightHandedOrientation = true;
else
    rightHandedOrientation = false;
end
