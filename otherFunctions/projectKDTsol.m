function [M, F, varCol] = projectKDTsol(varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.


%% Extract all needed data from options and varCol

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
noDofs = varCol.noCtrlPts;
extraGP = varCol.extraGP;

if varCol.solveForPtot
    analytic = @(x) varCol.analytic(x) + varCol.p_inc(x);
else
    analytic = varCol.analytic;
end

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
k = varCol.k;
totNoQP = 0;

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1);
[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP); 
Qxi = Q2D(:,1);
Qeta = Q2D(:,2);

%% Build global matrices
sizeMe = n_en^2;
spIdxRow = zeros(sizeMe,noElems);
spIdxCol = zeros(sizeMe,noElems);

Mvalues = zeros(sizeMe,noElems); 
Fvalues   = zeros(n_en,noElems,no_funcs); 
F_indices = zeros(n_en,noElems); 

%% Find nodes at which to evaluate the exact solution
%         for e = 1:noElems
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
    m_e = zeros(n_en);

    for gp_y = 1:size(W2D,1)
        pt_y = Q2D(gp_y,:);
        wt_y = W2D(gp_y);

        xi_y   = parent2ParametricSpace(Xi_e_y,  pt_y(1));
        eta_y  = parent2ParametricSpace(Eta_e_y, pt_y(2));

        [R_y, dRdxi_y, dRdeta_y] = NURBS2DBasis(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);
        J_y = pts_y'*[dRdxi_y' dRdeta_y'];
        crossProd_y = cross(J_y(:,1),J_y(:,2));
        J_1_y = norm(crossProd_y);
        n_y = crossProd_y/J_1_y;

        y = R_y*pts_y;
        fact_y = J_1_y * J_2_y * wt_y;
        y_d_n = P_far*n_y./norm2(P_far);
        x_d_y = P_far*y.'./norm2(P_far);
        dPhiFarField = -1i*k*y_d_n.*exp(-1i*k*x_d_y)/(4*pi);
        m_e = m_e + R(i,:)'*R(i,:) * abs(J_1(i)) * J_2 * W2D(i);  

        for e_z = 1:noElems   
            e_z
            [z,n_z,fact_z] = getBEMquadPtsKDT(e_y,e_z,W2D_2,Q,W,y,pt_y(1),pt_y(2),...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM);
            p_h_gp = 4*p_inc(z);
            p_h_gp(or(n_z*d_vec > 0,repmat(sum(n_z.*(y-z) > 0,2),1,size(d_vec,2)))) = 0;
            ymz = repmat(y,size(z,1),1) - z;
            r = norm2(ymz);

            if plotFarField
                p_h = p_h + (dPhiFarField * (dPhi_kdny(ymz,r,n_y.',k).* fact_z).'*p_h_gp).' * fact_y;  
            else
%                         xmy = P_far - repmat(x,size(P_far,1),1);
%                         r = norm2(xmy);
%                         p_h = p_h + p_h_gp.'.*dPhi_kdny(xmy,r,n_x.',k)* fact_x;  
            end
        end
    end
    spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
    spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
    Mvalues(:,e) = reshape(m_e, sizeMe, 1);

    F_indices(:,e) = sctr';
    Fvalues(:,e,:) = R.'*(p_h(:,:,e).*repmat(abs(J_1)*J_2.*W2D,1,no_funcs));
end        

%% Collect data into global matrices (and load vector)
F = vectorAssembly(Fvalues,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);


varCol.totNoQP = totNoQP;

