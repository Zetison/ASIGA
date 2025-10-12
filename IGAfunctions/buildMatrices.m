function task = buildMatrices(task,i_varCol)
% Create IGA global matrices
% Implemented for linear elasticity operator and the laplace operator, with
% possibility of computing the mass matrix and loading vector from body
% force. Thus, the function handles elasticity, laplace- and poisson equation, 
% and dynamic versions of these.
if nargin < 2
    i_varCol = 1;
end

%% Extract all needed data from options and varCol
degree = task.varCol{i_varCol}.degree;
knotVecs = task.varCol{i_varCol}.knotVecs;
index = task.varCol{i_varCol}.index;
noElems = task.varCol{i_varCol}.noElems;
elRange = task.varCol{i_varCol}.elRange;
element = task.varCol{i_varCol}.element;
element2 = task.varCol{i_varCol}.element2;
weights = task.varCol{i_varCol}.weights;
controlPts = task.varCol{i_varCol}.controlPts;
pIndex = task.varCol{i_varCol}.pIndex;
noDofs = task.varCol{i_varCol}.noDofs;
noCtrlPts = task.varCol{i_varCol}.noCtrlPts;

d_f = task.varCol{i_varCol}.fieldDimension;
d_p = task.varCol{i_varCol}.patches{1}.nurbs.d_p;
buildStiffnessMatrix = task.varCol{i_varCol}.buildStiffnessMatrix;
buildMassMatrix = task.varCol{i_varCol}.buildMassMatrix;
applyBodyLoading = task.varCol{i_varCol}.applyBodyLoading;
operator = task.varCol{i_varCol}.operator;

if strcmp(operator,'linearElasticity')
    C = task.varCol{i_varCol}.C;
else
    C = NaN; % Will not be used
end

extraGP = task.misc.extraGP;
isPML = task.varCol{i_varCol}.isPML;
usePML = any(isPML(:)) && isfield(task,'pml');
method = task.misc.method;

if usePML
    pml = task.pml;
    k = task.misc.omega/task.varCol{i_varCol}.c_f;
    if isnan(pml.gamma)
        switch pml.sigmaType
            case 1
                pml.gamma = 5;
            case 2
                pml.gamma = (pml.n+1)/(k*pml.t);
            case {3,4}
                pml.gamma = 1/(k*pml.t);
            case 5
                pml.gamma = pml.alpha/k/pml.n/pml.t.^(pml.n-1);
        end
    end
    formulation = task.misc.formulation;
    r_a = task.misc.r_a;
    task.pml = pml;
    d = task.pml.t;
else
    k = NaN;
    pml = NaN;
    r_a = NaN;
    d = NaN;
    isPML = false(noElems,1);
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
if usePML
    Kvalues = complex(Kvalues);
    Mvalues = complex(Mvalues);
end
if applyBodyLoading
    F_indices = zeros(d_f*n_en,noElems); 
    Fvalues = zeros(d_f*n_en,noElems); 
    force = task.varCol{i_varCol}.force;
else
    force = NaN;
end

[Q, W] = gaussTensorQuad(degree+1+extraGP(1:d_p));
if isfield(task.misc,'progressBars')
    progressBars = task.misc.progressBars;
else
    progressBars = false;
end
nProgressStepSize = ceil(noElems/1000);
if progressBars
    try
        ppm = ParforProgMon('Building mass/stiffness matrix: ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
else
    ppm = NaN;
end

%% Build global matrices
% for e = 1:noElems
parfor e = 1:noElems
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
    J = cell(1,3);
    if buildStiffnessMatrix
        if any(isPML(e,:)) && usePML && strcmp(formulation,'GSB')
            if pml.linearAbsorption
                switch sum(isPML(e,:))
                    case 1
                        j1 = find(isPML(e,:));
                        Isigma = intSigmaPML(xi(:,j1),pml);
                        for i = 1:d_p
                            J{i} = R{i+1}*pts;
                            if isPML(e,i)
                                J{i} = J{i}.*(1 + 1i*sigmaPML(xi(:,j1),pml));
                            else
                                d2X = (R{i+1}.*R{j1+1}./R{1})*pts;
                                J{i} = J{i} + 1i*d2X.*Isigma;
                            end
                        end
                    case 2
                        i_nonAbsorption = find(~isPML(e,:));
                        j_absorption = setdiff(1:d_p,i_nonAbsorption);
                        d2Xabsorption = (R{j_absorption(1)+1}.*R{j_absorption(2)+1}./R{1})*pts;
                        d3X = (R{2}.*R{3}.*R{4}./R{1}.^2)*pts;
                        Isigma = intSigmaPML(xi,pml);
                        sigma = sigmaPML(xi,pml);
                        for i = 1:d_p
                            dX = R{i+1}*pts;
                            J{i} = dX;

                            if i == i_nonAbsorption
                                for j2 = j_absorption
                                    d2X = (R{i+1}.*R{j2+1}./R{1})*pts;
                                    J{i} = J{i} + 1i*Isigma(:,j2).*d2X;
                                end
                                J{i} = J{i} - (1i*xi(:,j_absorption(1)).*Isigma(:,j_absorption(2)) + 1i*xi(:,j_absorption(2)).*Isigma(:,j_absorption(1)) - Isigma(:,j_absorption(1)).*Isigma(:,j_absorption(2))).*d3X;
                            else
                                J{i} = J{i} + 1i*sigma(:,i).*dX;
                                j1 = setdiff(j_absorption,i);
                                J{i} = J{i} - sigma(:,i).*(1i*xi(:,j1) - Isigma(:,j1)).*d2Xabsorption;
                            end
                        end
                    case 3
                        d3X = (R{2}.*R{3}.*R{4}./R{1}.^2)*pts;
                        Isigma = intSigmaPML(xi,pml);
                        sigma = sigmaPML(xi,pml);
                        xi1iI = xi + 1i*Isigma;
                        for i = 1:d_p
                            dX = R{i+1}*pts;
                            J{i} = dX;
                            otherIdx = setdiff(1:d_p,i);
                            temp = complex(zeros(size(J{i})));
                            for j2 = otherIdx
                                temp = temp - xi1iI(:,j2).*(R{i+1}.*R{j2+1}./R{1})*pts;
                            end
                            temp = temp + prod(xi1iI(:,otherIdx),2).*d3X;
                            J{i} = J{i} + 1i*sigma(:,i).*(dX + temp);
                        end
                end
            else
                xi_t = complex(xi);
    
                for i = 1:d_p
                    if isPML(e,i)
                        xi_t(:,i) = xi(:,i) + 1i*intSigmaPML(xi(:,i),pml);
                    end
                end
    
                R_t = NURBSbasis(I, xi_t, degree, knots, wgts);
                for i = 1:d_p
                    J{i} = R_t{i+1}*pts;
                    if isPML(e,i)
                        J{i} = J{i}.*(1+1i*sigmaPML(xi(:,i),pml));
                    end
                end
            end
        else
            for i = 1:d_p
                J{i} = R{i+1}*pts;
            end
        end
        J_1 = sum(J{1}.*cross(J{2},J{3},2),2);
    else
        J_1 = getJacobian(R,pts,d_p);
    end
    fact = J_1 * J_2.* W;
    useNIL = strcmp(method, 'NIL') && strcmp(operator, 'Laplace') && any(isPML(e,:));
    if useNIL
        X = R{1}*pts;
        x = X(:,1);
        y = X(:,2);
        z = X(:,3);
        rho = norm2(X);
        Om = 1 - (r_a-rho).^4/d^4;
        eL = 1 - (r_a-rho).^3.*(r_a+3*rho)/d^4;
        Omstr2 = -12*(r_a-rho).^2/d^4;
        
        A11 = (Om.^2.*x.^2 + eL.^2.*(y.^2+z.^2))./(eL.*rho.^2);
        A12 = ((Om.^2-eL.^2).*x.*y)./(eL.*rho.^2);
        A13 = ((Om.^2-eL.^2).*x.*z)./(eL.*rho.^2);
        A22 = (Om.^2.*y.^2 + eL.^2.*(x.^2+z.^2))./(eL.*rho.^2);
        A23 = ((Om.^2-eL.^2).*y.*z)./(eL.*rho.^2);
        A33 = (Om.^2.*z.^2 + eL.^2.*(x.^2+y.^2))./(eL.*rho.^2);
        B   = 2*1i*k*(eL - Om.^2)./eL;        
        C_NIL  = (Om.^2.*Omstr2.*rho.^2.*(1 - 1i*k*rho) + 2*eL.^2.*(1i + k*rho).*(1i*Om + k*rho) - eL.*Om.^2.*(-2 + k*rho.*(4i + k*rho)))./(eL.^2.*rho.^2);
    end

    if buildStiffnessMatrix        
        a11 = J{1}(:,1);
        a21 = J{1}(:,2);
        a31 = J{1}(:,3);
        a12 = J{2}(:,1);
        a22 = J{2}(:,2);
        a32 = J{2}(:,3);
        a13 = J{3}(:,1);
        a23 = J{3}(:,2);
        a33 = J{3}(:,3);
        JinvT1 = [(a22.*a33-a23.*a32)./J_1, (a23.*a31-a21.*a33)./J_1, (a21.*a32-a22.*a31)./J_1]; % First row of transpose of Jinv
        JinvT2 = [(a13.*a32-a12.*a33)./J_1, (a11.*a33-a13.*a31)./J_1, (a12.*a31-a11.*a32)./J_1]; % Second row of transpose of Jinv
        JinvT3 = [(a12.*a23-a13.*a22)./J_1, (a13.*a21-a11.*a23)./J_1, (a11.*a22-a12.*a21)./J_1]; % Third row of transpose of Jinv
        dRdX = cell(d_f,1);
        if usePML && strcmp(formulation,'STD')
            X = R{1}*pts;
            [r, theta, phi] = evaluateProlateCoords(X,0);
            if any(isPML(e,:))
                rs = (r-r_a)/pml.t;
                intSigma = pml.t*intSigmaPML(rs,pml);
                D = 1 + 1i*[sigmaPML(rs,pml), intSigma./r, intSigma./r];
            else
                D = ones(size(xi));
            end
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
        else
            dRdX{1} = JinvT1(:,1).*R{2} + JinvT1(:,2).*R{3} + JinvT1(:,3).*R{4};
            dRdX{2} = JinvT2(:,1).*R{2} + JinvT2(:,2).*R{3} + JinvT2(:,3).*R{4};
            dRdX{3} = JinvT3(:,1).*R{2} + JinvT3(:,2).*R{3} + JinvT3(:,3).*R{4};
        end

        if useNIL
            k_e = zeros(n_en);
            for i = 1:numel(fact)
                dRdx = [dRdX{1}(i,:); dRdX{2}(i,:); dRdX{3}(i,:)];
                Aij = [A11(i),A12(i),A13(i);
                       A12(i),A22(i),A23(i);
                       A13(i),A23(i),A33(i)];
                Bj = X(i,:)*B(i)/rho(i);
                k_e = k_e + (Aij*dRdx).'*dRdx * fact(i); 
                k_e = k_e - ((Bj*dRdx).'*R{1}(i,:)).' * fact(i); 
            end
            Kvalues(:,e) = reshape(k_e, (d_f*n_en)^2, 1);
        else
            Kvalues(:,e) = stiffnessElementMatrix(dRdX,fact,d_f,n_en,operator,C);
        end
    end
    if buildMassMatrix
        if useNIL
            Mvalues(:,e) = kron2(R{1},R{1}) * (fact.*C_NIL/k^2);
        else
            Mvalues(:,e) = kron2(R{1},R{1}) * fact;
        end
    end
    if applyBodyLoading
        v = R{1}*pts;
        if nargin(force) == 2
            [~,crossProd] = getJacobian(R,pts,d_p);
            n = crossProd./norm2(crossProd);
            f_gp = force(v,n);
        elseif nargin(force) == 3
            f_gp = force(v,e == noElems,J);
        else
            f_gp = force(v);
        end
        Fvalues(:,e) = kron2(f_gp,R{1}) * fact;
        F_indices(:,e) = sctr_k_e;
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
    task.varCol{i_varCol}.A_K = sparse(double(spIdxRow),double(spIdxCol),Kvalues,noDofs,noDofs,numel(Kvalues));
end

if applyBodyLoading
    task.varCol{i_varCol}.FF = vectorAssembly(Fvalues,F_indices,noDofs);
end

if buildMassMatrix
    task.varCol{i_varCol}.A_M = kron(sparse(double(spIdxRowM),double(spIdxColM),Mvalues,noCtrlPts,noCtrlPts,numel(Mvalues)),eye(d_f));
end

function sigma = sigmaPML(xi,pml)
sigma = zeros(size(xi));
switch pml.sigmaType
    case 0
        return
    case 1
        sigma = xi.*(exp(pml.gamma*xi)-1);
    case 2
        sigma = -pml.gamma*xi.^pml.n*log(pml.eps);
    case 3
        sigma = pml.gamma./(1-xi).^pml.n;
    case 4
        sigma = pml.gamma*(1./(1-xi).^pml.n - 1);
    case 5
        sigma = pml.gamma*xi.^pml.n;
    otherwise
        error('Not implemented')
end

function I = intSigmaPML(xi,pml)
% I = int_0^xi sigma(zeta) dzeta
I = zeros(size(xi));
switch pml.sigmaType
    case 0
        return
    case 1 % sigma(xi) = xi*(exp(gamma*xi)-1)
        gamma = pml.gamma;
        I = (exp(gamma*xi).*(gamma*xi-1) + 1)/gamma^2 - xi.^2/2;
    case 2 % sigma(xi) = gamma*xi^n
        n = pml.n;
        I = -pml.gamma*xi.^(n+1)/(n+1)*log(pml.eps);
    case 3 % sigma(xi) = gamma/(1-xi)^n
        n = pml.n;
        if n == 1
            I = -pml.gamma*log(1-xi);
        else
            I = pml.gamma*((1-xi).^(1-n)-1)./(n-1);
        end
    case 4 % sigma(xi) = gamma*(1/(1-xi)^n - 1)
        n = pml.n;
        if n == 1
            I = -pml.gamma*(log(1-xi) + xi);
        else
            I = pml.gamma*(((1-xi).^(1-n)-1)./(n-1) - xi);
        end
    case 5 % sigma(xi) = gamma*xi^n
        n = pml.n;
        I = pml.gamma*xi.^(n+1)/(n+1);
    otherwise
        error('Not implemented')
end






