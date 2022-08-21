function I = edgeLengths(nurbs,dirs,noExtraEvalPts)

if nargin < 3
    noExtraEvalPts = 1;
end
% This routine computes the length of all element edges and assumes no
% patches have internal C^-1 continuities
[Q, W] = gaussTensorQuad(12);
noPatches = numel(nurbs);
I = cell(1,noPatches);
for patch = 1:noPatches
    d_p = nurbs{patch}.d_p;
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
                noElems = numel(uniqueKnots)-1;
                for e = 1:noElems
                    Xi_e = uniqueKnots(e:e+1);
        
                    J_2 = (Xi_e(2)-Xi_e(1))/2;
                    xi_dir = parent2ParametricSpace(Xi_e, Q);
                    xi = repmat(parm_pts(i,:),numel(W),1);
                    xi(:,dir) = xi_dir;
                    X = cell(1,d_p+1);
                    [X{:}] = evaluateNURBS(nurbs{patch}, xi, 1);
                    I_xi = I_xi + norm2(X{dir+1}).' * J_2 * W;
                end
            end
            if I_xi > I_max
                I_max = I_xi;
            end
        end
        I{patch} = I_max;
    end
end









