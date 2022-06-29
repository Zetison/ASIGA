function I = edgeLengths(nurbs,dirs)
% This routine computes the length of all element edges and assumes no
% patches have internal C^-1 continuities
[Q, W] = gaussTensorQuad(12);
noPatches = numel(nurbs);
I = cell(1,noPatches);
for patch = 1:noPatches
    noElemsDir = zeros(1,3);
    for i = 1:d_p
        noElemsDir(i) = numel(unique(nurbs{patch}.knots{i}))-1;
    end
    I{patch} = zeros(noElemsDir);
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
    end
end









