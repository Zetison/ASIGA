function [M, F] = bestApproximation(varCol)
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
noDofs = varCol.noDofs;
extraGP = varCol.extraGP;

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
k = varCol.k;

switch varCol.formulation
    case 'SL2E'
        %% Preallocation and initiallizations
        n_en = (p_xi+1)*(p_eta+1);
        [W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP); 
        v_values   = zeros(size(W2D,1),noElems,3); 


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

            v_e = zeros(size(W2D,1),3);

            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));

                R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                v_e(gp,:) = R*pts;
            end
            v_values(:,e,:) = v_e;
        end
        v_values = reshape(v_values, size(W2D,1)*noElems, 3);
        analytic_values = varCol.analytic(v_values);
        analytic_values = reshape(analytic_values, size(W2D,1), noElems);


        %% Build global matrices
        sizeMe = n_en^2;
        spIdxRow = zeros(sizeMe,noElems);
        spIdxCol = zeros(sizeMe,noElems);

        Mvalues = zeros(sizeMe,noElems); 
        F_values   = zeros(n_en,noElems); 
        F_indices = zeros(n_en,noElems); 

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
            
            F_e = zeros(n_en,1);
            m_e = zeros(n_en);

            for gp = 1:size(W2D,1)
                pt = Q2D(gp,:);
                wt = W2D(gp);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));

                [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

                J = pts'*[dRdxi' dRdeta'];
                crossProd = cross(J(:,1),J(:,2));
                J_1 = norm(crossProd);

                y = R*pts;
                if useEnrichedBfuns
                    R = R*exp(1i*k*(y*d_vec));
                end
                m_e = m_e + R'*R * abs(J_1) * J_2 * wt;  

                p_analytic = analytic_values(gp,e);
                F_e = F_e + p_analytic*R' * abs(J_1) * J_2 * wt; 
            end
            spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
            spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
            Mvalues(:,e) = reshape(m_e, sizeMe, 1);

            F_indices(:,e) = sctr';
            F_values(:,e) = F_e;
        end
    case 'VL2E'
        elRangeZeta = varCol.elRange{3};
        p_zeta = varCol.degree(3); % assume p_eta is equal in all patches
        
        %% Preallocation and initiallizations
        n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);
        v_values = zeros(n_en,noElems,3); 

        [W3D,Q3D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP,p_zeta+1+extraGP); 

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

            v_e = zeros(size(W3D,1),3);

            for gp = 1:size(W3D,1)
                pt = Q3D(gp,:);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));
                zeta  = parent2ParametricSpace(Zeta_e, pt(3));

                R = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

                v_e(gp,:) = R*pts;
            end
            v_values(:,e,:) = v_e;
        end
        v_values = reshape(v_values, size(W3D,1)*noElems, 3);
        analytic_values = varCol.analytic(v_values);
        analytic_values = reshape(analytic_values, size(W3D,1), noElems);


        %% Build global matrices
        sizeMe = n_en^2;
        spIdxRow = zeros(sizeMe,noElems);
        spIdxCol = zeros(sizeMe,noElems);

        Mvalues = zeros(sizeMe,noElems); 
        F_values   = zeros(n_en,noElems); 
        F_indices = zeros(n_en,noElems); 

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
            
            J_2 = 0.125*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1))*(Zeta_e(2)-Zeta_e(1));

            F_e = zeros(n_en,1);
            m_e = zeros(n_en);

            for gp = 1:size(W3D,1)
                pt = Q3D(gp,:);
                wt = W3D(gp);

                xi   = parent2ParametricSpace(Xi_e,  pt(1));
                eta  = parent2ParametricSpace(Eta_e, pt(2));
                zeta  = parent2ParametricSpace(Zeta_e, pt(3));

                [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

                J = pts'*[dRdxi' dRdeta' dRdzeta'];
                J_1 = det(J);

                y = R*pts;
                if useEnrichedBfuns
                    R = R*exp(1i*k*dot(d_vec, y));
                end
                m_e = m_e + R.'*R * abs(J_1) * J_2 * wt;  

                p_analytic = analytic_values(gp,e);
                F_e = F_e + p_analytic*R.' * abs(J_1) * J_2 * wt; 
            end
            spIdxRow(:,e) = reshape(repmat(sctr',1,n_en),n_en^2,1);
            spIdxCol(:,e) = reshape(repmat(sctr,n_en,1),n_en^2,1);
            Mvalues(:,e) = reshape(m_e, sizeMe, 1);

            F_indices(:,e) = sctr';
            F_values(:,e) = F_e;
        end
end
        

%% Collect data into global matrices (and load vector)
F = vectorAssembly(F_values,F_indices,noDofs);

M = sparse(spIdxRow,spIdxCol,Mvalues,noDofs,noDofs);



