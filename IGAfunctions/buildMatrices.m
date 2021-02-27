function varCol = buildMatrices(varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles elasticity, laplace- and poisson equation, 
% and dynamic versions of these.

%% Extract all needed data from options and varCol
degree = varCol.degree;
knotVecs = varCol.knotVecs;
index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
pIndex = varCol.pIndex;
noDofs = varCol.noDofs;
noCtrlPts = varCol.noCtrlPts;

d_f = varCol.fieldDimension;
d_p = varCol.patches{1}.nurbs.d_p;
buildStiffnessMatrix = varCol.buildStiffnessMatrix;
buildMassMatrix = varCol.buildMassMatrix;
applyBodyLoading = varCol.applyBodyLoading;
operator = varCol.operator;

if strcmp(operator,'linearElasticity')
    C = varCol.C;
else
    C = NaN; % Will not be used
end

extraGP = varCol.extraGP;
Q = gaussTensorQuad(degree(1:2)+1+extraGP(1:d_p-1));

if isfield(varCol,'gamma')
    gamma = varCol.gamma;
    sigmaType = varCol.sigmaType;
    formulation = varCol.formulation;
    r_PML = varCol.r_PML;
    r_a = varCol.r_a;
    [zeta0Nodes_a, noElems_a, element_a, element2_a, index_a, pIndex_a] = meshBoundary(varCol,false);
    J_a = zeros(size(Q,1),3,3,noElems_a);
    for e = 1:noElems_a
    % parfor e = 1:noElems
        patch = pIndex_a(e);
        knots = knotVecs{patch}(1:2);
        Xi_e = zeros(2,2);
        for i = 1:2
            Xi_e(i,:) = elRange{i}(index_a(e,i),:);
        end

        sctr = zeta0Nodes_a(element_a(e,:));
        pts = controlPts(sctr,:);
        wgts = weights(element2_a(e,:),:);
        
        xi = parent2ParametricSpace(Xi_e, Q);
        I = findKnotSpans(degree(1:2), xi(1,:), knots);
        R = NURBSbasis(I, xi, degree(1:2), knots, wgts);
        J_a(:,:,1,e) = R{2}*pts;
        J_a(:,:,2,e) = R{3}*pts;
    end
%     J_a = permute(J_a,[1,3,2,4]);
else
    r_PML = NaN;
    r_a = NaN;
    gamma = 0;
    sigmaType = 0;
    decayDirs = zeros(noElems,3,'logical');
    formulation = 'GSB';
end

%% Preallocation and initiallizations
n_en = prod(degree+1);
sizeKe = (d_f*n_en)^2;
sizeMe = n_en^2;
spIdxRow = zeros(sizeKe,noElems,'uint32');
spIdxCol = zeros(sizeKe,noElems,'uint32');
spIdxRowM = zeros(sizeMe,noElems,'uint32');
spIdxColM = zeros(sizeMe,noElems,'uint32');
if buildStiffnessMatrix
    Kvalues = zeros(sizeKe,noElems); 
else
    Kvalues = NaN;
end

if buildMassMatrix
    Mvalues = zeros(sizeMe,noElems); 
else
    Mvalues = NaN;
end
if strcmp(varCol.method,'PML')
    Kvalues = complex(Kvalues);
    Mvalues = complex(Mvalues);
end
if applyBodyLoading
    F_indices = zeros(d_f*n_en,noElems); 
    Fvalues = zeros(d_f*n_en,noElems); 
    force = varCol.force;
else
    force = NaN;
end

[Q, W] = gaussTensorQuad(degree+1+extraGP(1:d_p));

progressBars = varCol.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    ppm = ParforProgMon('Building mass/stiffness matrix: ', noElems, nProgressStepSize);
else
    ppm = NaN;
end

%% Build global matrices
for e = 1:noElems
% parfor e = 1:noElems
	if progressBars && mod(e,nProgressStepSize) == 0
        ppm.increment();
	end
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(d_p,2);
    for i = 1:d_p
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;
    
    sctr_k_e = zeros(1,d_f*n_en);
    for i = 1:d_f
        sctr_k_e(i:d_f:end) = d_f*(sctr-1)+i;
    end
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts);
    J1 = R{2}*pts;
    J2 = R{3}*pts;
    J3 = R{4}*pts;
    if strcmp(formulation,'GSB')
        D = sigmaPML(xi(:,3),gamma,sigmaType);
        e_a = mod(e-1,noElems_a)+1;
        J1 = J1 + 1i*(J1 - repmat(J_a(:,:,1,e_a),degree(3)+1+extraGP(3),1)).*D(:,1);
        J2 = J2 + 1i*(J2 - repmat(J_a(:,:,2,e_a),degree(3)+1+extraGP(3),1)).*D(:,2);
        J3 = J3 + 1i*(J3 - repmat(J_a(:,:,3,e_a),degree(3)+1+extraGP(3),1)).*D(:,3);
        J_1 = sum(J1.*cross(J2,J3,2),2);
    else
        J_1 = getJacobian(R,pts,d_p);
    end
    fact = J_1 * J_2.* W;

    if buildStiffnessMatrix        
        a11 = J1(:,1);
        a21 = J1(:,2);
        a31 = J1(:,3);
        a12 = J2(:,1);
        a22 = J2(:,2);
        a32 = J2(:,3);
        a13 = J3(:,1);
        a23 = J3(:,2);
        a33 = J3(:,3);
        JinvT1 = [(a22.*a33-a23.*a32)./J_1, (a23.*a31-a21.*a33)./J_1, (a21.*a32-a22.*a31)./J_1]; % First row of transpose of Jinv
        JinvT2 = [(a13.*a32-a12.*a33)./J_1, (a11.*a33-a13.*a31)./J_1, (a12.*a31-a11.*a32)./J_1]; % Second row of transpose of Jinv
        JinvT3 = [(a12.*a23-a13.*a22)./J_1, (a13.*a21-a11.*a23)./J_1, (a11.*a22-a12.*a21)./J_1]; % Third row of transpose of Jinv
        dRdX = cell(d_f,1);
        switch formulation
            case 'GSB'
                dRdX{1} = JinvT1(:,1).*R{2} + JinvT1(:,2).*R{3} + JinvT1(:,3).*R{4};
                dRdX{2} = JinvT2(:,1).*R{2} + JinvT2(:,2).*R{3} + JinvT2(:,3).*R{4};
                dRdX{3} = JinvT3(:,1).*R{2} + JinvT3(:,2).*R{3} + JinvT3(:,3).*R{4};
            case 'STD'
                X = R{1}*pts;
                [r, theta, phi] = evaluateProlateCoords(X,0);
                rs = (r-r_PML)/(r_a-r_PML);
                intSigma = (r_a-r_PML)*intSigmaPML(rs,gamma,sigmaType);
                D = 1 + 1i*[rs.*exp(gamma*rs), intSigma./r, intSigma./r];
                sint = sin(theta);
                cost = cos(theta);
                sinp = sin(phi);
                cosp = cos(phi);
                Jsinv1 = [ sint.*cosp,      sint.*sinp,      cost];
                Jsinv2 = [ cost.*cosp./r,   cost.*sinp./r,  -sint./r];
                Jsinv3 = [-sinp./(r.*sint), cosp./(r.*sint), zeros(numel(W),1)];
                DinvJsinv1 = 1./D(:,1).*Jsinv1;
                DinvJsinv2 = 1./D(:,2).*Jsinv2;
                DinvJsinv3 = 1./D(:,3).*Jsinv3;
                JsDinvJsinv1 = sint.*cosp.*DinvJsinv1 + r.*cost.*cosp.*DinvJsinv2 - r.*sint.*sinp.*DinvJsinv3;
                JsDinvJsinv2 = sint.*sinp.*DinvJsinv1 + r.*cost.*sinp.*DinvJsinv2 + r.*sint.*cosp.*DinvJsinv3;
                JsDinvJsinv3 = cost.*DinvJsinv1       - r.*sint.*DinvJsinv2;
                JinvJsDinvJsJinv1 = JinvT1(:,1).*JsDinvJsinv1 + JinvT2(:,1).*JsDinvJsinv2 + JinvT3(:,1).*JsDinvJsinv3;
                JinvJsDinvJsJinv2 = JinvT1(:,2).*JsDinvJsinv1 + JinvT2(:,2).*JsDinvJsinv2 + JinvT3(:,2).*JsDinvJsinv3;
                JinvJsDinvJsJinv3 = JinvT1(:,3).*JsDinvJsinv1 + JinvT2(:,3).*JsDinvJsinv2 + JinvT3(:,3).*JsDinvJsinv3;
                dRdX{1} = JinvJsDinvJsJinv1(:,1).*R{2} + JinvJsDinvJsJinv2(:,1).*R{3} + JinvJsDinvJsJinv3(:,1).*R{4};
                dRdX{2} = JinvJsDinvJsJinv1(:,2).*R{2} + JinvJsDinvJsJinv2(:,2).*R{3} + JinvJsDinvJsJinv3(:,2).*R{4};
                dRdX{3} = JinvJsDinvJsJinv1(:,3).*R{2} + JinvJsDinvJsJinv2(:,3).*R{3} + JinvJsDinvJsJinv3(:,3).*R{4};
                fact =  fact.*prod(D,2);
        end

        Kvalues(:,e) = stiffnessElementMatrix(dRdX,fact,d_f,n_en,operator,C);
    end
    if buildMassMatrix
        Mvalues(:,e) = kron2(R{1},R{1}) * fact;
    end
    if applyBodyLoading
        v = R{1}*pts;
        f_gp = force(v);
        F_indices(:,e) = sctr_k_e;
        Fvalues(:,e) = kron2(f_gp,R{1}) * fact;
    end
    if buildStiffnessMatrix
        spIdxRow(:,e) = kron(ones(1,d_f*n_en),sctr_k_e);
        spIdxCol(:,e) = kron(sctr_k_e,ones(1,d_f*n_en));
    end
    if buildMassMatrix
        spIdxRowM(:,e) = kron(ones(1,n_en),sctr);
        spIdxColM(:,e) = kron(sctr,ones(1,n_en));
    end
end

%% Collect data into global matrices (and load vector)
if buildStiffnessMatrix
    varCol.A_K = sparse(double(spIdxRow),double(spIdxCol),Kvalues,noDofs,noDofs,numel(Kvalues));
end

if applyBodyLoading
    varCol.F = vectorAssembly(Fvalues,F_indices,noDofs);
end

if buildMassMatrix
    varCol.A_M = kron(sparse(double(spIdxRowM),double(spIdxColM),Mvalues,noCtrlPts,noCtrlPts,numel(Mvalues)),eye(d_f));
end

function sigma = sigmaPML(zeta,gamma,sigmaType)
sigma = zeros(numel(zeta),3);
switch sigmaType
    case 0
        return
    case 1
%         intSigma = s*intSigmaPML(zeta,gamma,sigmaType)./(1+s*zeta);
        intSigma = intSigmaPML(zeta,gamma,sigmaType)./zeta;
        sigma = [intSigma,intSigma,zeta.*exp(gamma*zeta)];
    otherwise
        error('Not implemented')
end

function I = intSigmaPML(zeta,gamma,sigmaType)
I = zeros(size(zeta));
switch sigmaType
    case 0
        return
    case 1
        I = (exp(gamma*zeta).*(gamma*zeta-1)+1)/gamma^2-zeta.^2/2;
    otherwise
        error('Not implemented')
end






