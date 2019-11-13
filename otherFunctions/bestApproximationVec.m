function [M, F, varCol] = bestApproximationVec(varCol)
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

switch varCol.formulation
    case 'SL2E'
        %% Preallocation and initiallizations
        n_en = (p_xi+1)*(p_eta+1);
        [W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP); 
        Qxi = Q2D(:,1);
        Qeta = Q2D(:,2);
        v_values   = zeros(size(W2D,1),noElems,3); 
        n_values   = zeros(size(W2D,1),noElems,3); 
        J_1_values   = zeros(size(W2D,1),noElems,1); 


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

            xi   = parent2ParametricSpace(Xi_e,  Qxi);
            eta  = parent2ParametricSpace(Eta_e, Qeta);
            [R, dRdxi, dRdeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
                
            J1 = dRdxi*pts;
            J2 = dRdeta*pts;
            crossProd = cross(J1,J2,2);
            J_1 = norm2(crossProd);
            n_values(:,e,:) = crossProd./J_1;
            v_values(:,e,:) = R*pts;
            J_1_values(:,e) = J_1;
        end
        v_values = reshape(v_values, size(W2D,1)*noElems, 3);
        n_values = reshape(n_values, size(W2D,1)*noElems, 3);
        if nargin(analytic) == 2
            analytic_values = analytic(v_values,n_values);
        else
            analytic_values = analytic(v_values);
        end
            
        no_funcs = size(analytic_values,2);
        analytic_values = permute(reshape(analytic_values, size(W2D,1), noElems, no_funcs),[1,3,2]);


        %% Build global matrices
        sizeMe = n_en^2;
        spIdxRow = zeros(sizeMe,noElems);
        spIdxCol = zeros(sizeMe,noElems);

        Mvalues = zeros(sizeMe,noElems); 
        Fvalues   = zeros(n_en,noElems,no_funcs); 
        F_indices = zeros(n_en,noElems); 

%         for e = 1:noElems
        parfor e = 1:noElems
            patch = pIndex(e); % New
            Xi = knotVecs{patch}{1}; % New
            Eta = knotVecs{patch}{2}; % New

            idXi = index(e,1);
            idEta = index(e,2);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);
            
            J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New
            
            
            xi   = parent2ParametricSpace(Xi_e,  Qxi);
            eta  = parent2ParametricSpace(Eta_e, Qeta);

            R = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

%             J1 = dRdxi*pts;
%             J2 = dRdeta*pts;
%             crossProd = cross(J1,J2,2);
%             J_1 = norm2(crossProd);
            J_1 = J_1_values(:,e);

            y = R*pts;
            if useEnrichedBfuns
                temp = exp(1i*k*(y*d_vec));
                R = R.*temp(:,ones(1,noGp));
            end
            m_e = zeros(n_en);
            for i = 1:numel(W2D)
                m_e = m_e + R(i,:)'*R(i,:) * abs(J_1(i)) * J_2 * W2D(i);  
            end
            
            spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
            spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
            Mvalues(:,e) = reshape(m_e, sizeMe, 1);

            F_indices(:,e) = sctr';
            Fvalues(:,e,:) = R.'*(analytic_values(:,:,e).*repmat(abs(J_1)*J_2.*W2D,1,no_funcs));
            totNoQP = totNoQP + numel(Qxi);
        end
    case 'VL2E'
        elRangeZeta = varCol.elRange{3};
        p_zeta = varCol.degree(3); % assume p_eta is equal in all patches
        
        %% Preallocation and initiallizations
        n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);

        [W3D,Q3D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP,p_zeta+1+extraGP); 
        Qxi = Q3D(:,1);
        Qeta = Q3D(:,2);
        Qzeta = Q3D(:,3);

        v_values = zeros(size(W3D,1),noElems,3); 
        %% Find nodes at which to evaluate the exact solution
%         for e = 1:noElems
        parfor e = 1:noElems
            patch = pIndex(e); % New
            Xi = knotVecs{patch}{1}; % New
            Eta = knotVecs{patch}{2}; % New
            Zeta = knotVecs{patch}{3}; % New

            idXi = index(e,1);
            idEta = index(e,2);
            idZeta = index(e,3);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);
            Zeta_e = elRangeZeta(idZeta,:);

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New

            xi   = parent2ParametricSpace(Xi_e,  Qxi);
            eta  = parent2ParametricSpace(Eta_e, Qeta);
            zeta = parent2ParametricSpace(Zeta_e, Qzeta);
            R = NURBS3DBasisVec(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);
                
            v_values(:,e,:) = R*pts;
        end
        v_values = reshape(v_values, size(W3D,1)*noElems, 3);
        analytic_values = analytic(v_values);
        no_funcs = size(analytic_values,2);
        analytic_values = permute(reshape(analytic_values, size(W3D,1), noElems, no_funcs),[1,3,2]);


        %% Build global matrices
        sizeMe = n_en^2;
        spIdxRow = zeros(sizeMe,noElems);
        spIdxCol = zeros(sizeMe,noElems);

        Mvalues = zeros(sizeMe,noElems); 
        Fvalues   = zeros(n_en,noElems,no_funcs); 
        F_indices = zeros(n_en,noElems); 

%         for e = 1:noElems
        parfor e = 1:noElems
            patch = pIndex(e); % New
            Xi = knotVecs{patch}{1}; % New
            Eta = knotVecs{patch}{2}; % New
            Zeta = knotVecs{patch}{3}; % New

            idXi = index(e,1);
            idEta = index(e,2);
            idZeta = index(e,3);

            Xi_e = elRangeXi(idXi,:);
            Eta_e = elRangeEta(idEta,:);
            Zeta_e = elRangeZeta(idZeta,:);

            J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

            sctr = element(e,:);
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:); % New


            xi   = parent2ParametricSpace(Xi_e,   Qxi);
            eta  = parent2ParametricSpace(Eta_e,  Qeta);
            zeta = parent2ParametricSpace(Zeta_e, Qzeta);
                
            [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVec(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

            dXdxi = dRdxi*pts;
            dXdeta = dRdeta*pts;
            dXdzeta = dRdzeta*pts;
            J_1 = dot(dXdxi,cross(dXdeta,dXdzeta,2),2);

            y = R*pts;
            if useEnrichedBfuns
                temp = exp(1i*k*(y*d_vec));
                R = R.*temp(:,ones(1,noGp));
            end
            m_e = zeros(n_en);
            for i = 1:numel(W3D)
                m_e = m_e + R(i,:)'*R(i,:) * abs(J_1(i)) * J_2 * W3D(i);  
            end

            
            spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
            spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
            Mvalues(:,e) = reshape(m_e, sizeMe, 1);

            F_indices(:,e) = sctr';
            Fvalues(:,e,:) = R.'*(analytic_values(:,:,e).*repmat(abs(J_1)*J_2.*W3D,1,no_funcs));
            totNoQP = totNoQP + numel(Qxi);
        end
end
        

%% Collect data into global matrices (and load vector)
F = vectorAssembly(Fvalues,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);


varCol.totNoQP = totNoQP;

