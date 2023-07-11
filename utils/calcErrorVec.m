function [relL2Error, relH1Error, relH1sError, relEnergyError] = calcErrorVec(task,i_o)

noDomains = length(task.varCol);
nodes = cell(noDomains,1);
factors = cell(noDomains,1);
u_hs = cell(noDomains,1);
du_hs = cell(noDomains,1);
for i = 1:noDomains
    if nargin < 1
        U = task.varCol{i}.U;
    else
        U = task.varCol{i}.U(:,i_o);
    end
    %% Extract all needed data from task.varCol
    if i == 1 && strcmp(task.misc.method,'PML')
        varCol_dummy = task.varCol{1};
        domain = varCol_dummy.geometry.domains.domain;
        [j,setFound] = findSet(domain,'Omega_a');
        if ~setFound
            error('Omega_a is not defined')
        end
        varCol_dummy.nurbs = varCol_dummy.nurbs(str2num(domain{j}.item{1}.Text));
        [j,setFound] = findSet(varCol_dummy.geometry.topologysets.set,'homDirichlet');
        if setFound
            varCol_dummy.geometry.topologysets.set(j) = [];
        end
        varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
        index = varCol_dummy.index;
        noElems = varCol_dummy.noElems;
        elRange = varCol_dummy.elRange;
        element = varCol_dummy.element;
        element2 = varCol_dummy.element2;
        weights = varCol_dummy.weights;
        controlPts = varCol_dummy.controlPts;
        knotVecs = varCol_dummy.knotVecs;
        pIndex = varCol_dummy.pIndex;
    else
        index = task.varCol{i}.index;
        noElems = task.varCol{i}.noElems;
        elRange = task.varCol{i}.elRange;
        element = task.varCol{i}.element;
        element2 = task.varCol{i}.element2;
        weights = task.varCol{i}.weights;
        controlPts = task.varCol{i}.controlPts;
        knotVecs = task.varCol{i}.knotVecs;
        pIndex = task.varCol{i}.pIndex;
    end
    d_p = task.varCol{i}.patches{1}.nurbs.d_p;

    degree = task.varCol{i}.degree; % assume degree is equal in all patches


    noCtrlPts = task.varCol{i}.noCtrlPts;

    %% Preallocation and initiallizations
    if false
        noQuadPts = 6;
        [W,Q] = gaussianQuadNURBS(noQuadPts,noQuadPts,noQuadPts); 
    else
        extraGP = task.misc.extraGP;
        [Q, W] = gaussTensorQuad(degree+3+extraGP);
    end
    if mod(i,2) == 0
        noComponents = 3;
        noComponentsDeriv = 9;
%         noComponentsDeriv = 6;
        
        Ux = U(1:3:3*noCtrlPts);
        Uy = U(2:3:3*noCtrlPts);
        Uz = U(3:3:3*noCtrlPts);
    else
        Ux = NaN;
        Uy = NaN;
        Uz = NaN;
        noComponents = 1;
        noComponentsDeriv = 3;
    end
    u_h = zeros(size(W,1),noElems,noComponents);
    du_h = zeros(size(W,1),noElems,noComponentsDeriv);
    fact = zeros(size(W,1),noElems);
    points = zeros(size(W,1),noElems,3);
    n_en = prod(degree+1);

%     for e = 1:noElems
    parfor e = 1:noElems
        patch = pIndex(e);
        knots = knotVecs{patch};
        Xi_e = zeros(d_p,2);
        for ii = 1:d_p
            Xi_e(ii,:) = elRange{ii}(index(e,ii),:);
        end

        sctr = element(e,:);
        if isnan(Ux)
            U_sctr = U(sctr,:);
        else
            Uxyz = [Ux(sctr) Uy(sctr) Uz(sctr)];
        end
        pts = controlPts(sctr,:);
        wgts = weights(element2(e,:),:);

        J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;

        xi = [parent2ParametricSpace(Xi_e, Q), zeros(size(Q,1),1)];
        I = findKnotSpans(degree, xi(1,:), knots);
        R = NURBSbasis(I, xi, degree, knots, wgts);
        J_1 = getJacobian(R,pts,d_p);

        dXdxi = R{2}*pts;
        dXdeta = R{3}*pts;
        dXdzeta = R{4}*pts;

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
        dRdx = repmat(Jinv1(:,1),1,n_en).*R{2} + repmat(Jinv1(:,2),1,n_en).*R{3} + repmat(Jinv1(:,3),1,n_en).*R{4};
        dRdy = repmat(Jinv2(:,1),1,n_en).*R{2} + repmat(Jinv2(:,2),1,n_en).*R{3} + repmat(Jinv2(:,3),1,n_en).*R{4};
        dRdz = repmat(Jinv3(:,1),1,n_en).*R{2} + repmat(Jinv3(:,2),1,n_en).*R{3} + repmat(Jinv3(:,3),1,n_en).*R{4};
        if mod(i,2) == 0
            u_h(:,e,:) = R{1}*Uxyz;

            dUdx = dRdx*Uxyz;
            dUdy = dRdy*Uxyz;
            dUdz = dRdz*Uxyz;

            du_h(:,e,:) = [dUdx, dUdy, dUdz];
        else
            u_h(:,e,:) = R{1}*U_sctr;
            du_h(:,e,:) = [dRdx*U_sctr,dRdy*U_sctr,dRdz*U_sctr];
        end
        fact(:,e) = J_1 * J_2 .* W;
        points(:,e,:) = R{1}*pts;
    end
    u_hs{i} = reshape(u_h, size(u_h,1)*size(u_h,2),noComponents);
    du_hs{i} = reshape(du_h, size(du_h,1)*size(du_h,2),noComponentsDeriv);
    factors{i} = reshape(fact, size(fact,1)*size(fact,2),1);
    nodes{i} = reshape(points, size(points,1)*size(points,2),3);
end

layer = task.analyticFunctions(nodes);

omega = task.misc.omega;

H1Error = 0;
H1sError = 0;
EnergyError = 0;
L2Error = 0;
normalizationEnergy = 0;
normalizationH1 = 0;
normalizationH1s = 0;
normalizationL2 = 0;
for i = 1:noDomains  
    switch layer{i}.media
        case 'solid'
            u = [layer{i}.u{1},layer{i}.u{2},layer{i}.u{3}];
            du = [layer{i}.du{1,1},layer{i}.du{1,2},layer{i}.du{1,3},...
                  layer{i}.du{2,1},layer{i}.du{2,2},layer{i}.du{2,3},...
                  layer{i}.du{3,1},layer{i}.du{3,2},layer{i}.du{3,3}];
            sigma = [layer{i}.sigma{1},layer{i}.sigma{2},layer{i}.sigma{3},layer{i}.sigma{4},layer{i}.sigma{5},layer{i}.sigma{6}];

            C = task.varCol{i}.C;
            strain_vec = (C\sigma.').';
            strain_h = [du_hs{i}(:,1),du_hs{i}(:,5),du_hs{i}(:,9),du_hs{i}(:,6)+du_hs{i}(:,8),du_hs{i}(:,7)+du_hs{i}(:,3),du_hs{i}(:,2)+du_hs{i}(:,4)];
            strain_e = strain_vec-strain_h;

            eCe = real(sum((strain_e*C).*conj(strain_e),2)); % the usage of real() is to remove machine epsilon imaginary part
            uCu = real(sum((strain_vec*C).*conj(strain_vec),2)); % the usage of real() is to remove machine epsilon imaginary part

            u_e = u-u_hs{i};
            du_e = du-du_hs{i};

            u2 = sum(abs(u).^2,2);
            du2 = sum(abs(du).^2,2);
            u_e2 = sum(abs(u_e).^2,2);
            du_e2 = sum(abs(du_e).^2,2);

            H1Error = H1Error + sum((u_e2 + du_e2).*factors{i});
            H1sError = H1sError + sum(du_e2.*factors{i});
            normalizationH1 = normalizationH1 + sum((u2 + du2).*factors{i});
            normalizationH1s = normalizationH1s + sum(du2.*factors{i});

            EnergyError = EnergyError + sum((eCe + task.varCol{i}.rho*omega^2*u_e2).*factors{i});
            normalizationEnergy = normalizationEnergy + sum((uCu + task.varCol{i}.rho*omega^2*u2).*factors{i});

            L2Error = L2Error + sum(u_e2.*factors{i});
            normalizationL2 = normalizationL2 + sum(u2.*factors{i});
        case 'fluid'
            p = layer{i}.p;
            dp = [layer{i}.dp{1}, layer{i}.dp{2}, layer{i}.dp{3}];
            p_e = p-u_hs{i};
            dp_e = dp-du_hs{i};

            p2 = abs(p).^2;
            dp2 = sum(abs(dp).^2,2);
            p_e2 = abs(p_e).^2;
            dp_e2 = sum(abs(dp_e).^2,2);

            H1Error = H1Error + sum((p_e2 + dp_e2).*factors{i});
            H1sError = H1sError + sum(dp_e2.*factors{i});
            normalizationH1 = normalizationH1 + sum((p2 + dp2).*factors{i});
            normalizationH1s = normalizationH1s + sum(dp2.*factors{i});
            k = task.misc.omega/task.varCol{i}.c_f;

            EnergyError = EnergyError + 1/(task.varCol{i}.rho*omega^2)*sum((dp_e2 + k^2*p_e2).*factors{i});
            normalizationEnergy = normalizationEnergy + 1/(task.varCol{i}.rho*omega^2)*sum((dp2 + k^2*p2).*factors{i});

            L2Error = L2Error + sum(p_e2.*factors{i});
            normalizationL2 = normalizationL2 + sum(p2.*factors{i});
    end
end

relEnergyError = 100*sqrt(EnergyError/normalizationEnergy);
relH1Error = 100*sqrt(H1Error/normalizationH1);
relH1sError = 100*sqrt(H1sError/normalizationH1s);
relL2Error = 100*sqrt(L2Error/normalizationL2);