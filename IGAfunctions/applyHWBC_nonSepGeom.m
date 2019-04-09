function F = applyHWBC_nonSepGeom(varCol,alpha_s_arr)

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

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

dp_inc = varCol.dp_inc;

extraGP = varCol.extraGP;

noDofs_tot = varCol.noDofs_tot;


% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.
noDofs = varCol.noDofs;
D = varCol.D;
N = varCol.N;
F = zeros(noDofs_tot,length(alpha_s_arr));        % external force vector

n_en = (p_xi+1)*(p_eta+1);
Fvalues = zeros(n_en,noElems,length(alpha_s_arr));
indices = zeros(n_en,noElems);
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP); 
parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:)); % New   
    

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    f_e = zeros(n_en,length(alpha_s_arr));
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R_fun, dRxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

        J = pts'*[dRxi' dRdeta'];
        crossProd = cross(J(:,1), J(:,2)); 

        v = R_fun*pts;
        
        n = -crossProd.'/norm(crossProd);
        
        deriv = -dp_inc(v,n).';
        
        f_e = f_e + R_fun'*deriv*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctr';
    Fvalues(:,e,:) = f_e;
end
% If IEbasisfuns satisfy Kronecker delta property on surface, we can do the
% following
% for alpha_s_Nr = 1:length(alpha_s_arr)
%     F(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
% end

for alpha_s_Nr = 1:length(alpha_s_arr)
    for i = 1:N       
        if sum(D(i,:)) ~= 0
            F(:,alpha_s_Nr) = F(:,alpha_s_Nr) + vectorAssembly(sum(D(i,:))*Fvalues(:,:,alpha_s_Nr),indices+(i-1)*noDofs,noDofs_tot);
        end
    end
end
