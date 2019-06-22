function relError = calcSurfErrorBndryMethodVec(varCol, U, LpOrder)

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

% if noElems < 600
%     noQuadPts = 16;
% else
%     noQuadPts = 10;
% end
% noQuadPts = max(9,max(p_xi,p_eta)+3);
% noQuadPts = max(9,max(p_xi,p_eta)+3);
% [W2D,Q2D] = gaussianQuadNURBS(noQuadPts,noQuadPts);  
extraGP = varCol.extraGP;
% [W2D,Q2D] = gaussianQuadNURBS(p_xi+3+extraGP,p_eta+3+extraGP);  
[W2D,Q2D] = gaussianQuadNURBS(p_xi+10,p_eta+10);  

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

    xi   = parent2ParametricSpace(Xi_e,  Q2D(:,1));
    eta  = parent2ParametricSpace(Eta_e, Q2D(:,2));

    [R, dRdxi, dRdeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
    
    J1 = dRdxi*pts;
    J2 = dRdeta*pts;
    crossProd = cross(J1,J2,2);
    J_1 = norm2(crossProd);
        
    y = R*pts;
    if useEnrichedBfuns
        R = R*repmat(exp(1i*k*(y*d_vec)),1,size(R,2));
    end
    if ~isMFS
        U_sctr = U(sctr);
        p_h(:,e) = R*U_sctr;
    end
    fact(:,e) = J_1 * J_2 .* W2D;
    points(:,e,:) = y;
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

if varCol.solveForPtot
    analytic = @(x) varCol.analytic(x) + varCol.p_inc(x);
else
    analytic = varCol.analytic;
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


