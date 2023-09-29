function I = NURBSarcLength(nurbs,a,b,parm_pt,dir,pp1qp,noQP)
if nargin < 6
    pp1qp = false;
end
if nargin < 7
    noQP = 12;
end
if pp1qp
    [Q, W] = gaussTensorQuad(nurbs.degree(dir)+3);
else
    [Q, W] = gaussTensorQuad(noQP);
end
uniqueKnots = unique(nurbs.knots{dir});
I = 0;
d_p = nurbs.d_p;
for e = 1:numel(uniqueKnots)-1
    Xi_e = uniqueKnots(e:e+1);
    if Xi_e(2) <= a || b <= Xi_e(1)
        continue
    end
    if b < Xi_e(2)
        Xi_e(2) = b;
    end
    if Xi_e(1) < a
        Xi_e(1) = a;
    end
    J_2 = (Xi_e(2)-Xi_e(1))/2;
    xi_dir = parent2ParametricSpace(Xi_e, Q);
    xi = repmat(parm_pt,numel(W),1);
    xi(:,dir) = xi_dir;
    X = cell(1,d_p+1);
    [X{:}] = evaluateNURBS(nurbs, xi, 1);
    I = I + norm2(X{dir+1}).' * J_2 * W;
end
