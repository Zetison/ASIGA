function relError = calcH1norm_volume(varCol, U)

%% Extract all needed data from varCol
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);
p_zeta = varCol.nurbs.degree(3);

index = varCol.index;
noElems = varCol.noElems;
element = varCol.element;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
weights = varCol.weights;
controlPts = varCol.controlPts;


analytic = varCol.analytic;
gAnalytic = varCol.gAnalytic;

%% Preallocation and initiallizations
if noElems < 600
    noQuadPts = 16;
else
    noQuadPts = 10;
end
[W3D,Q3D] = gaussianQuadNURBS(noQuadPts,noQuadPts,noQuadPts); 

p_h = zeros(size(W3D,1),noElems);
gradient_h = zeros(size(W3D,1),noElems,3);
fact = zeros(size(W3D,1),noElems);
points = zeros(size(W3D,1),noElems,3);

parfor e = 1:noElems
% for e = 1:noElems
    idXi = index(e,1);
    idEta = index(e,2);
    idZeta = index(e,3);
    
    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);
    Zeta_e = elRangeZeta(idZeta,:);
    
    J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));
    
    p_h_temp = zeros(size(W3D,1),1);
    gradient_h_temp = zeros(size(W3D,1),3);
    fact_temp = zeros(size(W3D,1),1);
    points_temp = zeros(size(W3D,1),3);
    sctr = element(e,:);
    U_sctr = U(sctr);
    pts = controlPts(sctr,:);
    
    for gp = 1:size(W3D,1)
        pt = Q3D(gp,:);
        wt = W3D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));
        zeta = parent2ParametricSpace(Zeta_e,pt(3));
        
        [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);

        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        J_1 = det(J);
        X = R*pts; 
        p_h_temp(gp) = R*U_sctr;
        gradient_h_temp(gp,:) = J'\[dRdxi; dRdeta; dRdzeta]*U_sctr;
        
        fact_temp(gp) = J_1 * J_2 * wt;
        points_temp(gp,:) = X;
    end
    p_h(:,e) = p_h_temp;
    gradient_h(:,e,:) = gradient_h_temp;
    fact(:,e) = fact_temp;
    points(:,e,:) = points_temp;
end
p_h = reshape(p_h, size(p_h,1)*size(p_h,2),1);
gradient_h = reshape(gradient_h, size(gradient_h,1)*size(gradient_h,2),3);
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
p_exact = analytic(points);
gradient_exact = gAnalytic(points);

Error = sqrt(sum((abs(p_exact - p_h).^2 + norm2(gradient_exact-gradient_h).^2).*fact));
normalization = sqrt(sum((abs(p_exact).^2 + norm2(gradient_exact).^2).*fact));

relError = 100*Error/normalization;


