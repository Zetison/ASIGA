function I = edgeLengths(nurbs,dirs)
% This routine computes the length of all element edges and assumes no
% patches have internal C^-1 continuities
[Q, W] = gaussTensorQuad(12);
for patch = 1:numel(nurbs)
    for dir = dirs
        uniqueXi = unique(nurbs{patch}.knots{dir});
        dir2 = setdiff(1:d_p,dir);
        
        knots = nurbs{patch}.knots(dir2);
        for j = 1:numel(dir2)
            knots{j} = insertUniform(unique(knots{j}), noExtraEvalPts);
        end
        if d_p == 2
            uniqueEta = knots{1};
        elseif d_p == 3
            uniqueEta = [kron(knots{1},ones(numel(knots{2}),1)), kron(ones(numel(knots{1}),1),knots{2})];
        end
        parm_pts = NaN(size(uniqueEta,1),d_p);
        parm_pts(:,dir2) = uniqueEta;


        uniqueKnots = unique(nurbs{patch}.knots{dir});
for j = 1:numel(uniqueXi)-1
    I_max = -Inf;
    for i = 1:size(parm_pts,1)
        I = 0;
        d_p = nurbs{patch}.d_p;
        noElems = numel(uniqueKnots)-1;
        for e = 1:noElems
            Xi_e = uniqueKnots(e:e+1);

            J_2 = (Xi_e(2)-Xi_e(1))/2;
            xi_dir = parent2ParametricSpace(Xi_e, Q);
            xi = repmat(parm_pt,numel(W),1);
            xi(:,dir) = xi_dir;
            X = cell(1,d_p+1);
            [X{:}] = evaluateNURBS(nurbs{patch}, xi, 1);
            I = I + norm2(X{dir+1}).' * J_2 * W;
        end
    end
end
I = zeros(noElems,2*d_p);
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
    X = R{1}*pts;
    I(e,:) = sum(norm2(X{dir+1}).' * J_2 * W);
end











