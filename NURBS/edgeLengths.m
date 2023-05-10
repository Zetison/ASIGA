function edgeLen = edgeLengths(nurbs,dirs,noExtraEvalPts)
% This routine computes the maximum distance between two knots for all
% knots and for all directions in dirs for all patches in nurbs

if nargin < 3
    noExtraEvalPts = 2*max(nurbs{1}.degree)-1;
end
[Q, W] = gaussTensorQuad(12);
noQP = numel(W);
noPatches = numel(nurbs);
edgeLen = cell(1,noPatches);
for patch = 1:noPatches
    d_p = nurbs{patch}.d_p;
    d = nurbs{patch}.d;
    edgeLen{patch} = cell(1,d_p);
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
        noEtaEvalPts = size(uniqueEta,1);
        parm_pts = NaN(noEtaEvalPts,d_p);
        parm_pts(:,dir2) = uniqueEta;


        noElems = numel(uniqueXi)-1;
        XI = zeros(noQP*noEtaEvalPts,noElems,d_p);
        J_2 = (uniqueXi(2:end)-uniqueXi(1:end-1))/2;
        for e = 1:noElems
            Xi_e = uniqueXi(e:e+1);
            xi_dir = parent2ParametricSpace(Xi_e, Q);
            XI(:,e,:) = kron(parm_pts,ones(numel(W),1));
            XI(:,e,dir) = repmat(xi_dir,noEtaEvalPts,1);
        end
        X = cell(1,d_p+1);
        XI = reshape(XI,noQP*noEtaEvalPts*noElems,d_p);
        [X{:}] = evaluateNURBSvec(nurbs{patch}, XI, 1);
        X{dir+1} = reshape(X{dir+1},noQP,noEtaEvalPts*noElems,d);
        J_2 = repmat(kron(J_2,ones(1,noEtaEvalPts)),noQP,1);
        edgeLen{patch}{dir} = max(reshape(sum(vecnorm(X{dir+1},2,3) .* J_2 .* W, 1), noEtaEvalPts, noElems),[], 1);
    end
end








