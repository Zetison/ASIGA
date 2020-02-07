function [K, M, F] = buildGlobalMatricesVec(varCol, newOptions)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles both static and dynamic linear
% elasticity, laplace- and poisson equation, and dynamic versions of these.

%% Interpret input arguments
% set default options values
options = struct('operator','Laplace',...
                 'fieldDimension',1,...
                 'buildMassMatrix',0,...
                 'applyBodyLoading',0);
options = updateOptions(options,newOptions);

%% Extract all needed data from options and varCol
d = options.fieldDimension;
buildMassMatrix = options.buildMassMatrix;
applyBodyLoading = options.applyBodyLoading;
operator = options.operator;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
p_zeta = varCol.degree(3); % assume p_zeta is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
elRangeZeta = varCol.elRange{3};
pIndex = varCol.pIndex;
noDofs = varCol.noDofs;

if strcmp(operator,'linearElasticity')
    C = varCol.C;
else
    C = 0; % Will not be used
end

%% Preallocation and initiallizations
n_en = (p_xi+1)*(p_eta+1)*(p_zeta+1);
sizeKe = (d*n_en)^2;
spIdxRow = zeros(sizeKe,noElems,'uint32');
spIdxCol = zeros(sizeKe,noElems,'uint32');
Kvalues = zeros(sizeKe,noElems); 

if buildMassMatrix
    Mvalues = zeros(sizeKe,noElems); 
end
if applyBodyLoading
    F_indices = zeros(d*n_en,noElems); 
    Fvalues = zeros(d*n_en,noElems); 
    force = varCol.force;
else
    force = NaN;
end

extraGP = varCol.extraGP;
[W3D,Q3D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP,p_zeta+1+extraGP); 
Qxi = Q3D(:,1);
Qeta = Q3D(:,2);
Qzeta = Q3D(:,3);
% [Qxi, Wxi] = gaussQuad(p_xi+1);
% [Qeta, Weta] = gaussQuad(p_eta+1);
% [Qzeta, Wzeta] = gaussQuad(p_zeta+1);
% W3D = tensorWeights(Wxi,Weta,Wzeta);
%% Build global matrices
% warning('Not running in Parallel')
% keyboard
for e = 1:noElems
% parfor e = 1:noElems
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
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    k_e = zeros(d*n_en);
    if buildMassMatrix
        m_e = zeros(n_en);
    end
    if applyBodyLoading
        f_e = zeros(d*n_en,1);
    end

    xi   = parent2ParametricSpace(Xi_e,  Qxi);
    eta  = parent2ParametricSpace(Eta_e, Qeta);
    zeta = parent2ParametricSpace(Zeta_e,Qzeta);
    
%     [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVecTen(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);
    [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVec(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, wgts);

    dXdxi = dRdxi*pts;
    dXdeta = dRdeta*pts;
    dXdzeta = dRdzeta*pts;
    J_1 = dot(dXdxi,cross(dXdeta,dXdzeta,2),2);
                
    a11 = dXdxi(:,1);
    a21 = dXdxi(:,2);
    a31 = dXdxi(:,3);
    a12 = dXdeta(:,1);
    a22 = dXdeta(:,2);
    a32 = dXdeta(:,3);
    a13 = dXdzeta(:,1);
    a23 = dXdzeta(:,2);
    a33 = dXdzeta(:,3);
    Jinv1 = [(a22.*a33-a23.*a32)./J_1, (a23.*a31-a21.*a33)./J_1, (a21.*a32-a22.*a31)./J_1];
    Jinv2 = [(a13.*a32-a12.*a33)./J_1, (a11.*a33-a13.*a31)./J_1, (a12.*a31-a11.*a32)./J_1];
    Jinv3 = [(a12.*a23-a13.*a22)./J_1, (a13.*a21-a11.*a23)./J_1, (a11.*a22-a12.*a21)./J_1];
    dRdx = repmat(Jinv1(:,1),1,n_en).*dRdxi + repmat(Jinv1(:,2),1,n_en).*dRdeta + repmat(Jinv1(:,3),1,n_en).*dRdzeta;
    dRdy = repmat(Jinv2(:,1),1,n_en).*dRdxi + repmat(Jinv2(:,2),1,n_en).*dRdeta + repmat(Jinv2(:,3),1,n_en).*dRdzeta;
    dRdz = repmat(Jinv3(:,1),1,n_en).*dRdxi + repmat(Jinv3(:,2),1,n_en).*dRdeta + repmat(Jinv3(:,3),1,n_en).*dRdzeta;

    for i = 1:numel(W3D)
        dRdX = [dRdx(i,:); dRdy(i,:); dRdz(i,:)];
        switch operator
            case 'linearElasticity'
                B = strainDispMatrix3d(n_en,dRdX);
                k_e = k_e + B' * C * B * abs(J_1(i)) * J_2 * W3D(i); 
            case 'Laplace'
                k_e = k_e + dRdX.'*dRdX * abs(J_1(i)) * J_2 * W3D(i);  
        end
        if buildMassMatrix
            m_e = m_e + R(i,:)'*R(i,:) * abs(J_1(i)) * J_2 * W3D(i);  
        end

        if applyBodyLoading
            v = R(i,:)*pts;
            f_gp = force(v(1),v(2),v(3));
            f_e = f_e + kron(R(i,:).', f_gp) * abs(J_1(i)) * J_2 * W3D(i);
        end
    end

%     spIdxRow(:,e) = copyVector(sctr_k_e,d*n_en,1);
%     spIdxCol(:,e) = copyVector(sctr_k_e,d*n_en,2);
    spIdxRow(:,e) = kron(ones(1,d*n_en),sctr_k_e);
    spIdxCol(:,e) = kron(sctr_k_e,ones(1,d*n_en));
    temp = zeros(d*n_en,d*n_en);
    for i = 1:d
        for j = 1:d
            temp(i:d:end, j:d:end) = k_e(1+(i-1)*n_en:i*n_en, 1+(j-1)*n_en:j*n_en);
        end
    end
    Kvalues(:,e) = reshape(temp, sizeKe, 1);
    
    if buildMassMatrix
        temp = zeros(d*n_en,d*n_en);
        for i = 1:d
            temp(i:d:end, i:d:end) = m_e;
        end
        Mvalues(:,e) = reshape(temp, sizeKe, 1);
    end
    if applyBodyLoading
        F_indices(:,e) = sctr_k_e;
        Fvalues(:,e) = f_e;
    end
end

%% Collect data into global matrices (and load vector)
if applyBodyLoading
    F = vectorAssembly(Fvalues,F_indices,noDofs);
end

optimization = 3; % at the expence of speed, optimization=1 reduces most memory consumption, optimization=2 reduces moderately memory consumption, optimization=3 does not reduce memory consumption

switch optimization
    case 1
        primes_noElems = factor(noElems);
        prt = prod(primes_noElems(1:3));
        Kvalues = reshape(Kvalues, sizeKe*prt,noElems/prt);
        Mvalues = reshape(Mvalues, sizeKe*prt,noElems/prt);
        spIdxRow = reshape(spIdxRow, sizeKe*prt,noElems/prt);
        spIdxCol = reshape(spIdxCol, sizeKe*prt,noElems/prt);
        for i = 1:noElems/prt
            [spIdx,~,I] = unique([spIdxRow(:,i), spIdxCol(:,i)],'rows','stable');
            KvaluesTemp = accumarray(I,Kvalues(:,i));
            if buildMassMatrix
                MvaluesTemp = accumarray(I,Mvalues(:,i));
            end
            m = numel(KvaluesTemp);
            if m < sizeKe*prt
                Kvalues(:,i) = [KvaluesTemp; zeros(sizeKe*prt-m,1)];
                spIdxRow(:,i) = [spIdx(:,1); zeros(sizeKe*prt-m,1)];
                spIdxCol(:,i) = [spIdx(:,2); zeros(sizeKe*prt-m,1)];
                if buildMassMatrix
                    Mvalues(:,i) = [MvaluesTemp; zeros(sizeKe*prt-m,1)];
                end
            end
        end
    %     if noElems > 260000
    %         keyboard
    %     end
        nnzEntries = sum(~(~spIdxRow(:)));
        spIdxRowUnique = zeros(nnzEntries,1,'uint32');
        spIdxColUnique = zeros(nnzEntries,1,'uint32');
        KvaluesUnique = zeros(nnzEntries,1);
        if buildMassMatrix
            MvaluesUnique = zeros(nnzEntries,1);
        end
        counter = 1;
        for i = 1:noElems/prt
            indices = find(spIdxRow(:,i));
            m = numel(indices);
            spIdxRowUnique(counter:counter+m-1) = spIdxRow(indices,i);
            spIdxColUnique(counter:counter+m-1) = spIdxCol(indices,i);
            KvaluesUnique(counter:counter+m-1) = Kvalues(indices,i);
            if buildMassMatrix
                MvaluesUnique(counter:counter+m-1) = Mvalues(indices,i);
            end
            counter = counter + m;
        end
        Kvalues = KvaluesUnique;
        clear spIdxRow spIdxCol KvaluesUnique
        if buildMassMatrix
            Mvalues = MvaluesUnique;
            clear MvaluesUnique
        end
        spIdx = [spIdxRowUnique, spIdxColUnique];
        clear spIdxRowUnique spIdxColUnique

        [spIdx,~,I] = unique(spIdx,'rows','stable');
        Kvalues = accumarray(I,Kvalues);
        if buildMassMatrix
            Mvalues = accumarray(I,Mvalues);
        end
    case 2
        spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
        spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
        Kvalues = reshape(Kvalues,numel(Kvalues),1);
        if buildMassMatrix
            Mvalues = reshape(Mvalues,numel(Mvalues),1);
        end
        spIdx = [spIdxRow, spIdxCol];
        clear spIdxRow spIdxCol
        [spIdx,~,I] = unique(spIdx,'rows','stable');
        Kvalues = accumarray(I,Kvalues);
        if buildMassMatrix
            Mvalues = accumarray(I,Mvalues);
        end
end

if optimization == 3
    K = sparse(double(spIdxRow),double(spIdxCol),Kvalues,noDofs,noDofs,numel(Kvalues));
    clear Kvalues

    if buildMassMatrix
        M = sparse(double(spIdxRow),double(spIdxCol),Mvalues,noDofs,noDofs,numel(Mvalues));
    else
        M = [];
    end
else
    K = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Kvalues,noDofs,noDofs,numel(Kvalues));
    clear Kvalues

    if buildMassMatrix
        M = sparse(double(spIdx(:,1)),double(spIdx(:,2)),Mvalues,noDofs,noDofs,numel(Mvalues));
    else
        M = [];
    end
end


