function relError = calcSurfErrorVec(varCol, U, LpOrder)

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
elRangeZeta = varCol.elRange{3};
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;

analytic = varCol.analytic;

% Find elements on the inner surface for evaluation of
% backscattered pressure in far field
surfaceElements = [];
for e = 1:noElems
    idZeta = index(e,3);
    Zeta_e = elRangeZeta(idZeta,:); % [zeta_k,zeta_k+1]                    
    if Zeta_e(1) == 0
        surfaceElements = [surfaceElements e];
    end
end

% if noElems < 600
%     noQuadPts = 16;
% else
%     noQuadPts = 10;
% end
[W2D,Q2D] = gaussianQuadNURBS(p_xi+3,p_eta+3);  

p_h = zeros(size(W2D,1),length(surfaceElements));
fact = zeros(size(W2D,1),length(surfaceElements));
points = zeros(size(W2D,1), length(surfaceElements),3);

% for i = 1:length(surfaceElements) %[8 7 3 4 5 6 2 1]% 
parfor i = 1:length(surfaceElements) %[8 7 3 4 5 6 2 1]% 

    e = surfaceElements(i);
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    Zeta = knotVecs{patch}{3}; % New
    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:); % [xi_i,xi_i+1]
    Eta_e = elRangeEta(idEta,:); % [eta_j,eta_j+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctr = element(e,:);
    U_sctr = U(sctr,:);

    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:); % New
        
    xi   = parent2ParametricSpace(Xi_e,  Q2D(:,1));
    eta  = parent2ParametricSpace(Eta_e, Q2D(:,2));
    zeta = zeros(size(W2D));

    [R, dRdxi, dRdeta] = NURBS3DBasisVec(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

    J1 = dRdxi*pts;
    J2 = dRdeta*pts;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
    
    p_h(:,i) = R*U_sctr;
    fact(:,i) = J_1 * J_2 .* W2D;
    points(:,i,:) = R*pts;
end
p_h = reshape(p_h, size(p_h,1)*size(p_h,2),1);
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
p = analytic(points);

if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*fact))^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*fact))^(1/LpOrder);
end

relError = 100*Error/normalization;


