function relError = calcSurfErrorBndryMethodVec(varCol, U, LpOrder)

i = 1;
p_xi = varCol{i}.degree(1); % assume p_xi is equal in all patches
p_eta = varCol{i}.degree(2); % assume p_eta is equal in all patches

index = varCol{i}.index;
noElems = varCol{i}.noElems;
elRangeXi = varCol{i}.elRange{1};
elRangeEta = varCol{i}.elRange{2};
element = varCol{i}.element;
element2 = varCol{i}.element2;
weights = varCol{i}.weights;
controlPts = varCol{i}.controlPts;
knotVecs = varCol{i}.knotVecs;
pIndex = varCol{i}.pIndex;

k = varCol{i}.k;

if strcmp(varCol{1}.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol{1}.d_vec;
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
extraGP = varCol{1}.extraGP;
[W2D,Q2D] = gaussianQuadNURBS(p_xi+3+extraGP,p_eta+3+extraGP);  
% [W2D,Q2D] = gaussianQuadNURBS(p_xi+10,p_eta+10);  

p_h = zeros(size(W2D,1), noElems);
fact = zeros(size(W2D,1), noElems);
points = zeros(size(W2D,1), noElems,3);
isMFS = strcmp(varCol{1}.method,'MFS');

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
        U_sctr = U{i}(sctr);
        p_h(:,e) = R*U_sctr;
    end
    fact(:,e) = J_1 * J_2 .* W2D;
    points(:,e,:) = y;
end
fact = reshape(fact, size(fact,1)*size(fact,2),1);
points = reshape(points, size(points,1)*size(points,2),3);
if isMFS
    Phi_k = @(r) exp(1i*k*r)./(4*pi*r);
    n_cp = numel(U{i});
    y_s = varCol{i}.y_s;
    r = zeros(size(points));
    for j = 1:n_cp
        r(:,j) = norm2(elementAddition(-y_s(j,:), points));
    end
    p_h = Phi_k(r)*U{i};
else
    p_h = reshape(p_h, size(p_h,1)*size(p_h,2),1);
end

if varCol{i}.solveForPtot && varCol{i}.exteriorProblem
    analytic = @(x) varCol{i}.analytic(x) + varCol{i}.p_inc(x);
else
    analytic = varCol{i}.analytic;
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


