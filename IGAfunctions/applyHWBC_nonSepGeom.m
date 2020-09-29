function F = applyHWBC_nonSepGeom(varCol,alpha_s_arr)

degree = varCol.degree; % assume degree is equal in all patches
index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
dp_inc = varCol.dp_inc;
extraGP = varCol.extraGP;
noDofs_tot = varCol.noDofs_tot;
N = varCol.N;
F = zeros(noDofs_tot,length(alpha_s_arr));        % external force vector

n_en = prod(degree+1);
Fvalues = zeros(n_en,noElems,length(alpha_s_arr));
indices = zeros(n_en,noElems);
[Q, W] = gaussTensorQuad(degree+1+extraGP);
parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(2,2);
    for i = 1:2
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:)); % New  
            
    R = NURBSbasis(I, xi, degree, knots, wgts);
    [J_1, crossProd] = getJacobian(R,pts,2);
    X = R{1}*pts;
    n = -crossProd./repmat(J_1,1,3);
    deriv = -dp_inc(X,n);
    indices(:,e) = sctr';    
    Fvalues(:,e,:) = R{1}'*(deriv.*J_1 * J_2 .* W);
end
if strcmp(varCol.IEbasis,'Standard')
    noDofs = varCol.noDofs;
    D = generateCoeffMatrix(varCol);

    for alpha_s_Nr = 1:length(alpha_s_arr)
        for i = 1:N       
            if sum(D(i,:)) ~= 0
                F(:,alpha_s_Nr) = F(:,alpha_s_Nr) + vectorAssembly(sum(D(i,:))*Fvalues(:,:,alpha_s_Nr),indices+(i-1)*noDofs,noDofs_tot);
            end
        end
    end
else
    for alpha_s_Nr = 1:length(alpha_s_arr)
        F(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
    end
end