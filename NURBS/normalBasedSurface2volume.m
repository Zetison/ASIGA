function [nurbsVol,nurbs2] = normalBasedSurface2volume(nurbs,t_pml)
if ~iscell(nurbs)
    nurbs = {nurbs};
end
if 1
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
    task.varCol{1}.force = @(X,n) X + t_pml*n;
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
        FF(end+1) = X(end,:).*normal + t_pml;

        U = A\FF;
        nurbsVol = loftNURBS({nurbs,nurbs1,nurbs2},2);
    end
else
    task.misc.method = 'BA';
    task.misc.progressBars = false;
    task.misc.extraGP = [0,0,0];
    task.varCol{1}.dimension = 3;
    task.varCol{1}.nurbs = nurbs(1);
    task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
    task.varCol{1}.operator = 'laplace';
    task.varCol{1}.applyBodyLoading = true;
    task.varCol{1}.buildMassMatrix = true;
    task.varCol{1}.buildStiffnessMatrix = false;
    task.varCol{1}.fieldDimension = 3;
    task.varCol{1}.force = @(X,n) X + t_pml*n;
    task.varCol{1}.extraGP = 3*ones(1,2);
    task.varCol{1}.progressBars = false;
    task.varCol{1}.media = 'fluid';
    
    %% Extract all needed data from options and varCol
    degree = task.varCol{1}.degree;
    knotVecs = task.varCol{1}.knotVecs;
    index = task.varCol{1}.index;
    noElems = task.varCol{1}.noElems;
    elRange = task.varCol{1}.elRange;
    element = task.varCol{1}.element;
    element2 = task.varCol{1}.element2;
    weights = task.varCol{1}.weights;
    controlPts = task.varCol{1}.controlPts;
    pIndex = task.varCol{1}.pIndex;
    noDofs = task.varCol{1}.noDofs;

    d_f = task.varCol{1}.fieldDimension;
    d_p = task.varCol{1}.patches{1}.nurbs.d_p;

    extraGP = task.misc.extraGP;
    nen = prod(degree+1);

    %% Preallocation and initiallizations
    n_en = prod(degree+1);
    sizeKe = (d_f*n_en)^2;
    spIdxRow = zeros(sizeKe,noElems,'uint32');
    spIdxCol = zeros(sizeKe,noElems,'uint32');
    Kvalues = zeros(sizeKe,noElems); 
    F_indices = zeros(d_f*n_en,noElems); 
    Fvalues = zeros(d_f*n_en,noElems); 

    [Q, W] = gaussTensorQuad(degree+1+extraGP(1:d_p));
    
    %% Build global matrices
    for e = 1:noElems
%     parfor e = 1:noElems
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
        
        [J_1,n] = getJacobian(R,pts,d_p);
        fact = J_1 * J_2.* W;
        X = R{1}*pts;

        k_e = zeros(d_f*n_en);
        temp = zeros(d_f*n_en);
        for i = 1:numel(fact)
            RR = R{1}(i,:).'*R{1}(i,:);
            temp(1:nen,1:nen)                 =  0;
            temp(1:nen,nen+1:2*nen)           =  n(i,3)*RR;
            temp(1:nen,2*nen+1:3*nen)         = -n(i,2)*RR;
            
            temp(nen+1:2*nen,1:nen)           = -n(i,3)*RR;
            temp(nen+1:2*nen,nen+1:2*nen)     =  0;
            temp(nen+1:2*nen,2*nen+1:3*nen)   =  n(i,1)*RR;
            
            temp(2*nen+1:3*nen,1:nen)         =  n(i,2)*RR;
            temp(2*nen+1:3*nen,nen+1:2*nen)   = -n(i,1)*RR;
            temp(2*nen+1:3*nen,2*nen+1:3*nen) =  0;
            
            k_e = k_e + temp * fact(i); 
        end
        for i = 1:d_f
            for j = 1:d_f
                temp(i:d_f:end, j:d_f:end) = k_e(1+(i-1)*n_en:i*n_en, 1+(j-1)*n_en:j*n_en);
            end
        end        
        Kvalues(:,e) = reshape(temp, (d_f*n_en)^2, 1);
        Fvalues(:,e) = kron2(cross(n,X),R{1}) * fact;
        F_indices(:,e) = sctr_k_e;

        spIdxRow(:,e) = kron(ones(1,d_f*n_en),sctr_k_e);
        spIdxCol(:,e) = kron(sctr_k_e,ones(1,d_f*n_en));
    end

    %% Collect data into global matrices (and load vector)
    task.varCol{1}.A_M = sparse(double(spIdxRow),double(spIdxCol),Kvalues,noDofs,noDofs,numel(Kvalues));
    task.varCol{1}.FF = vectorAssembly(Fvalues,F_indices,noDofs);
    
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
%     nurbs2 = cell(size(nurbs));
%     for patch = 1:numel(nurbs)
%         P = nurbs{patch}.coeffs(1:3,:,:);
%         number = nurbs{patch}.number;
%         Pt = zeros(size(P));
%         
%         % Compute values at corner points
%         n = cross(P(:,2,1)-P(:,1,1),P(:,1,2)-P(:,1,1));
%         Pt(:,1,1) = P(:,1,1) + t_pml*n/norm(n);
%         n = cross(P(:,end,2)-P(:,end,1),P(:,end-1,1)-P(:,end,1));
%         Pt(:,end,1) = P(:,end,1) + t_pml*n/norm(n);
%         n = cross(P(:,1,end-1)-P(:,1,end),P(:,2,end)-P(:,1,end));
%         Pt(:,1,end) = P(:,1,end) + t_pml*n/norm(n);
%         n = cross(P(:,end-1,end)-P(:,end,end),P(:,end,end-1)-P(:,end,end));
%         Pt(:,end,end) = P(:,end,end) + t_pml*n/norm(n);
%         
%         % Compute values at boundary
%         % south boundary
%         for i = 2:number(1)-1
%             n_avg = (cross(P(:,i+1,1)-P(:,i,1), P(:,i,2)-P(:,i,1)) + ...
%                      cross(P(:,i,2)-P(:,i,1), P(:,i-1,1)-P(:,i,1)))/2;
%             n = cross(P(:,i,2)-P(:,i,1), P(:,i-1,1)-P(:,i,1));
%             Pt(:,i,1) = P(:,i,1) + dot(Pt(:,i-1,1)-P(:,i,1), n)/dot(n_avg,n)*n_avg; 
%         end
%         % west boundary
%         for j = 2:number(2)-1
%             n_avg = (cross(P(:,2,j)-P(:,1,j), P(:,1,j+1)-P(:,1,j)) + ...
%                      cross(P(:,1,j-1)-P(:,1,j), P(:,2,j)-P(:,1,j)))/2;
%             n = cross(P(:,1,j-1)-P(:,1,j), P(:,2,j)-P(:,1,j));
%             Pt(:,1,j) = P(:,1,j) + dot(Pt(:,1,j-1)-P(:,1,j), n)/dot(n_avg,n)*n_avg; 
%         end
%         % north boundary
%         for i = 2:number(1)-1
%             n_avg = (cross(P(:,i+1,end)-P(:,i,end), P(:,i,end-1)-P(:,i,end)) + ...
%                      cross(P(:,i,end-1)-P(:,i,end), P(:,i-1,end)-P(:,i,end)))/2;
%             n = cross(P(:,i,end-1)-P(:,i,end), P(:,i-1,end)-P(:,i,end));
%             Pt(:,i,end) = P(:,i,end) + dot(Pt(:,i-1,end)-P(:,i,end), n)/dot(n_avg,n)*n_avg; 
%         end
%         % east boundary
%         for j = 2:number(2)-1
%             n_avg = (cross(P(:,end,j+1)-P(:,end,j), P(:,end-1,j)-P(:,end,j)) + ...
%                      cross(P(:,end-1,j)-P(:,end,j),P(:,end,j-1)-P(:,end,j)))/2;
%             n = cross(P(:,end-1,j)-P(:,end,j), P(:,end,j-1)-P(:,end,j));
%             Pt(:,end,j) = P(:,end,j) + dot(Pt(:,end,j-1)-P(:,end,j), n)/dot(n_avg,n)*n_avg; 
%         end
%         
%         % Compute interior points
%         for j = 2:number(2)-1
%             for i = 2:number(1)-1
%                 n_avg = (cross(P(:,i+1,j)-P(:,i,j), P(:,i,j+1)-P(:,i,j)) + ...
%                          cross(P(:,i,j+1)-P(:,i,j), P(:,i-1,j)-P(:,i,j)) + ...
%                          cross(P(:,i-1,j)-P(:,i,j), P(:,i,j-1)-P(:,i,j)) + ...
%                          cross(P(:,i,j-1)-P(:,i,j), P(:,i+1,j)-P(:,i,j)))/4;
%                 n = cross(P(:,i-1,j)-P(:,i,j), P(:,i,j-1)-P(:,i,j));
%                 Pt(:,i,j) = P(:,i,j) + dot(Pt(:,i-1,j)-P(:,i,j), n)/dot(n_avg,n)*n_avg;  
%             end
%         end
%         coeffs = nurbs{patch}.coeffs;
%         coeffs(1:3,:,:) = Pt;
%         nurbs2(patch) = createNURBSobject(coeffs,nurbs{patch}.knots);
%     end
%     nurbsVol = loftNURBS({nurbs,nurbs2});
end
                
                
        
        
        
        
        
        
        


