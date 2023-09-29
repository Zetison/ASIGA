function [isEqual,orient] = NURBSisEqual(nurbs1,nurbs2,withOrient)
if nargin < 3
    withOrient = false;
end
orient = 0;
isEqual = nurbs1.d_p == nurbs2.d_p && ...
          nurbs1.d == nurbs2.d && ...
          all(sort(nurbs1.number) == sort(nurbs2.number)) && ...
          all(sort(nurbs1.degree) == sort(nurbs2.degree));
if ~isEqual
    return
end
Eps = 1e-10;
d_p = nurbs1.d_p;
if withOrient
    [flips,indices] = getOrientPerms(d_p);
    for i = 1:size(indices,1)
        for j = 1:2^d_p
            coeffs2 = nurbs2.coeffs;
            knots2 = nurbs2.knots;
            for jj = 1:numel(flips{j})
                flipIdx = flips{j}(jj);
                coeffs2 = flip(coeffs2,flipIdx+1);
                knots2{flipIdx} = 1-flip(knots2{flipIdx});
            end
            if d_p > 0
                coeffs2 = permute(coeffs2,[1,indices(i,:)+1]);
            end
            knots2 = knots2(indices(i,:));
            equalKnotVecs = true;
            for jj = 1:d_p
                if any(nurbs1.number ~= nurbs2.number(indices(i,:))) || ...
                   any(nurbs1.degree ~= nurbs2.degree(indices(i,:))) || ...
                   any(abs(nurbs1.knots{jj} - knots2{jj}) > Eps)
                    equalKnotVecs = false;
                    break
                end
            end
            isEqual = all(abs(nurbs1.coeffs(:) - coeffs2(:)) < Eps) && equalKnotVecs;
            if isEqual
                return
            end
            orient = orient + 1;
        end
    end
else
    for i = 1:d_p
        if any(abs(nurbs1.knots{i} - nurbs2.knots{i}) > Eps)
            isEqual = false;
        end
    end
    if ~isEqual
        return
    end
    temp = abs(nurbs1.coeffs - nurbs2.coeffs) < Eps;
    isEqual = all(temp(:));
end