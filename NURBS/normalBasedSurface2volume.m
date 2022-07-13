function [nurbsVol,nurbs2] = normalBasedSurface2volume(nurbs,t_pml,X_bApprox,colMethod)
if ~iscell(nurbs)
    nurbs = {nurbs};
end
if nargin < 3
    X_bApprox = 'BA';
%     X_bApprox = 'interp';
end
if nargin < 4
    colMethod = 'Grev';
%     colMethod = 'CG';
%     colMethod = 'GL';
end

noPatches = numel(nurbs);
nurbs2 = cell(1,noPatches);
task.misc.method = 'BA';
task.misc.progressBars = false;
task.misc.extraGP = [0,0,0];
task.varCol{1}.nurbs = nurbs;
task.varCol{1}.fieldDimension = 3;
switch X_bApprox
    case 'BA'
        task.varCol{1}.dimension = 3;
        task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
        task.varCol{1}.operator = 'laplace';
        task.varCol{1}.applyBodyLoading = true;
        task.varCol{1}.buildMassMatrix = true;
        task.varCol{1}.buildStiffnessMatrix = false;
        task.varCol{1}.force = @(X,n) X + t_pml*n;
        task.varCol{1}.extraGP = 3*ones(1,2);
        task.varCol{1}.progressBars = false;
        task = buildMatrices(task,1);
    case 'interp'
        task.varCol{1}.dimension = 1;
        task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
        degree = task.varCol{1}.degree; % assume degree is equal in all patches
        
        element = task.varCol{1}.element;
        element2 = task.varCol{1}.element2;
        weights = task.varCol{1}.weights;
        controlPts = task.varCol{1}.controlPts;
        knotVecs = task.varCol{1}.knotVecs;
        noElemsPatch = task.varCol{1}.noElemsPatch;
        noPatches = task.varCol{1}.noPatches;
        dofsToRemove = task.varCol{1}.dofsToRemove;
        noDofs = task.varCol{1}.noDofs;
        patches = task.varCol{1}.patches;
        d_f = task.varCol{1}.fieldDimension;
        noCtrlPts = task.varCol{1}.noCtrlPts;

        [cp_p, patchIdx, n_cp] = getCollocationPts(patches,noDofs,dofsToRemove,1e9*eps,colMethod);
        n_en = prod(degree+1);
        values = zeros(n_cp,n_en);
        FF = zeros(n_cp,d_f);

        spIdxRowM = zeros(n_cp,n_en,'uint32');
        spIdxColM = zeros(n_cp,n_en,'uint32');
        parfor i = 1:n_cp
            patch = patchIdx(i);
            knots = knotVecs{patch};
            Xi_x = knotVecs{patch}{1}; % New
            Eta_x = knotVecs{patch}{2}; % New
            uniqueXi = unique(Xi_x);
            uniqueEta = unique(Eta_x);
            noElementsXi = length(uniqueXi)-1;
            noElementsEta = length(uniqueEta)-1;

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
            x = R{1}*pts;
            J_temp = [R{2}; R{3}]*pts;
            m_1 = J_temp(1,:);
            m_2 = J_temp(2,:);
            crossProd = cross(m_1,m_2);

            n = crossProd/norm(crossProd);


            values(i,:) = R{1};
            FF(i,:) = x + t_pml*n;
    
            spIdxRowM(i,:) = i;
            spIdxColM(i,:) = sctr;
        end

        task.varCol{1}.FF = FF;
    
        task.varCol{1}.A_M = sparse(double(spIdxRowM),double(spIdxColM),values,n_cp,noCtrlPts,numel(values));
        task.misc.method = 'interp';
        task.misc.formulation = 'C';
end
task.varCol{1}.media = 'fluid';
[task,FF,~,~,A2] = collectMatrices(task);
task.UU = A2\FF;
task = postProcessSolution(task);
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
                
                
        
        
        
        
        
        
        


