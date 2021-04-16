function [nurbsVol,nurbs2] = normalBasedSurface2volume(nurbs,D)
if ~iscell(nurbs)
    nurbs = {nurbs};
end
task.misc.method = 'BA';
task.misc.progressBars = false;
task.misc.extraGP = [0,0,0];
task.varCol{1}.dimension = 3;
task.varCol{1}.nurbs = nurbs;
task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
task.varCol{1}.operator = 'laplace';
task.varCol{1}.applyBodyLoading = true;
task.varCol{1}.buildMassMatrix = true;
task.varCol{1}.buildStiffnessMatrix = false;
task.varCol{1}.fieldDimension = 3;
task.varCol{1}.force = @(X,n) X + D*n;
task.varCol{1}.extraGP = 3*ones(1,2);
task.varCol{1}.progressBars = false;
task = buildMatrices(task,1);
task.varCol{1}.media = 'fluid';
[task,FF,~,~,A2] = collectMatrices(task);
UU = A2\FF;
task = postProcessSolution(task,UU);
noPatches = numel(nurbs);
nurbs2 = cell(1,noPatches);
counter = 1;
for i = 1:noPatches
    sz = size(nurbs{i}.coeffs);
    coeffs = zeros(sz);
    coeffs(4,:,:) = nurbs{i}.coeffs(4,:,:);
    noDofs = task.varCol{1}.noCtrlPtsPatch(i)*task.varCol{1}.dimension;
    coeffs(1:3,:,:) = reshape(task.varCol{1}.U(counter:counter+noDofs-1,:).',3,sz(2),sz(3));
    counter = counter + noDofs;
    nurbs2(i) = createNURBSobject(coeffs,nurbs{i}.knots);
end
nurbsVol = loftNURBS({nurbs,nurbs2});

if 0
    task.varCol{1}.dimension = 1;
    task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
    patches = task.varCol{1}.patches;
    element = task.varCol{1}.element;
    element2 = task.varCol{1}.element2;
    knotVecs = task.varCol{1}.knotVecs;
    dofsToRemove = task.varCol{1}.dofsToRemove;
    noDofs = task.varCol{1}.noDofs;
    degree = task.varCol{1}.degree;
    noElemsPatch = task.varCol{1}.noElemsPatch;
    controlPts = task.varCol{1}.controlPts;
    weights = task.varCol{1}.weights;
    
    n_cp = noDofs - length(dofsToRemove);
    counter2 = 1;
    counter = 1;
    cp_p = zeros(n_cp,2);
    patchIdx = zeros(n_cp,1);
    p_zeta = 2;

    for patch = 1:noPatches
        nurbs = patches{patch}.nurbs;
        n_xi = nurbs.number(1);
        n_eta = nurbs.number(2);
        Xi_y = nurbs.knots{1};
        Eta_y = nurbs.knots{2};

        for j = 1:n_eta
            eta_bar = sum(Eta_y(j+1:j+degree(2)))/degree(2);
            for i = 1:n_xi
                if ~any(dofsToRemove == counter)
                    xi_bar = sum(Xi_y(i+1:i+degree(1)))/degree(1);
                    cp_p(counter2,:) = [xi_bar, eta_bar];
                    patchIdx(counter2) = patch;
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    end
    n_en = prod(degree+1);
    A = zeros(n_cp, noDofs);
    FF = zeros(n_cp, 3);
    crossProds = zeros(n_cp, 3);
    X = zeros(n_cp, 3);
    for i = 1:n_cp
%     parfor i = 1:n_cp
        patch = patchIdx(i);
        knots = knotVecs{patch};
        Xi_x = knotVecs{patch}{1}; % New
        Eta_x = knotVecs{patch}{2}; % New
        uniqueXi = unique(Xi_x);
        uniqueEta = unique(Eta_x);
        noElementsXi = length(uniqueXi)-1;
        noElementsEta = length(uniqueEta)-1;

        A_row = zeros(1, noDofs);
        xi_x = cp_p(i,1);
        eta_x = cp_p(i,2);

        xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
        eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
        e = sum(noElemsPatch(1:patch-1)) + xi_idx + noElementsXi*(eta_idx-1);

        sctr = element(e,:);
        pts = controlPts(sctr,:);
        wgts = weights(element2(e,:),:); % New

        xi = [xi_x, eta_x];
        I = findKnotSpans(degree, xi(1,:), knots);
        R = NURBSbasis(I, xi, degree, knots, wgts);
        R_x = R{1};
        dR_xdxi = R{2};
        dR_xdeta = R{3};
        x = R_x*pts;
        J_temp = [dR_xdxi; dR_xdeta]*pts;
        m_1 = J_temp(1,:);
        m_2 = J_temp(2,:);
        crossProd = cross(m_1,m_2);
        for j = 1:n_en
            A_row(sctr(j)) = A_row(sctr(j)) + R_x(j);
        end     
        A(i,:) = A_row;
        FF(i,:) = x;
        crossProds(i,:) = crossProd;
        X(i,:) = x;
    end
    A = kron(A,eye(3));
    A(:,end+1) = reshape(crossProds.'/p_zeta,[],1);
    normal = crossProds(end,:)/norm(crossProds(end,:));
    A(end+1,1:3) = A(end,1:3).*normal;
    FF = reshape(FF,[],1);
    FF(end+1) = X(end,:).*normal + D;
    
    U = A\FF;
    nurbsVol = loftNURBS({nurbs,nurbs1,nurbs2},2);
end


