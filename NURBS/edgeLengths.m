function I = edgeLengths(nurbs,dirs)
% This routine computes the length of all element edges and assumes no
% patches have internal C^-1 continuities
[Q, W] = gaussTensorQuad(12);
for patch = 1:numel(nurbs)
    for dir = dirs
        uniqueKnots = unique(nurbs{patch}.knots{dir});
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
