function relError = calcSurfErrorBndryMethod(varCol, U, LpOrder)
error('Depricated use calcSurfErrorBndryMethodVec instead')

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;

k = varCol.k;

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end

analytic = varCol.analytic;

if noElems < 600
    noQuadPts = 16;
else
    noQuadPts = 10;
end
noQuadPts = max(9,max(p_xi,p_eta)+3);
[W2D,Q2D] = gaussianQuadNURBS(noQuadPts,noQuadPts);  

p_h = zeros(size(W2D,1), noElems);
fact = zeros(size(W2D,1), noElems);
points = zeros(size(W2D,1), noElems,3);
isMFS = strcmp(varCol.method,'MFS');

parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:); % New
    
    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    if ~isMFS
        p_h_temp = zeros(size(W2D,1),1);
        U_sctr = U(sctr);
    end
    fact_temp = zeros(size(W2D,1),1);
    points_temp = zeros(size(W2D,1),3);
    normals_temp = zeros(size(W2D,1),3);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi   = parent2ParametricSpace(Xi_e,  pt(1));
        eta  = parent2ParametricSpace(Eta_e, pt(2));

        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

        J = pts'*[dRdxi' dRdeta'];
        crossProd = cross(J(:,1),J(:,2));
        J_1 = norm(crossProd);
        n = -crossProd/J_1;
        
        y = R*pts;
        if useEnrichedBfuns
            R = R*exp(1i*k*(y*d_vec));
        end
        if ~isMFS
            p_h_temp(gp) = R*U_sctr;
        end

        fact_temp(gp) = norm(crossProd) * J_2 * wt;
        
        points_temp(gp,:) = y;
        normals_temp(gp,:) = n;

    end
    if ~isMFS
        p_h(:,e) = p_h_temp;
    end
    fact(:,e) = fact_temp;
    points(:,e,:) = points_temp;
end
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
if isMFS
    Phi_k = @(r) exp(1i*k*r)./(4*pi*r);
    n_cp = numel(U);
    y_s = varCol.y_s;
    r = zeros(size(points));
    for j = 1:n_cp
        r(:,j) = norm2(elementAddition(-y_s(j,:), points));
    end
    p_h = Phi_k(r)*U;
else
    p_h = reshape(p_h, size(p_h,1)*size(p_h,2),1);
end

if strcmp(varCol.method, 'BEM') && ~strcmp(varCol.applyLoad, 'pointPulsation')
    p_inc = varCol.p_inc;
    p_h = p_h - p_inc(points);
end
p = analytic(points);
if isinf(LpOrder)
    Error = max(abs(p - p_h));
    normalization = max(abs(p));
else
    Error = (sum((abs(p - p_h).^LpOrder).*fact))^(1/LpOrder);
    normalization = (sum((abs(p).^LpOrder).*fact))^(1/LpOrder);
end
relError = 100*Error/normalization;


