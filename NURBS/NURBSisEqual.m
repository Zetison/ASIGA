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
    flips = cell(2^d_p,1);
    flips{1} = [];      % For d_p >= 1
    flips{2} = 1;       % For d_p >= 1
    flips{3} = 2;       % For d_p >= 2
    flips{4} = [1,2];   % For d_p >= 2
    flips{5} = 3;       % For d_p >= 3
    flips{6} = [1,3];   % For d_p >= 3
    flips{7} = [2,3];   % For d_p >= 3
    flips{8} = [1,2,3]; % For d_p >= 3
    indices = perms(1:d_p);
    for i = 1:size(indices,1)
        for j = 1:2^d_p
            coeffs2 = nurbs2.coeffs;
            knots = nurbs2.knots;
            for jj = 1:numel(flips{j})
                flipIdx = flips{j}(jj);
                coeffs2 = flip(coeffs2,flipIdx+1);
                knots{flipIdx} = 1-flip(knots{flipIdx});
            end
            temp = coeffs2;
            if d_p > 0
                temp = permute(coeffs2,[1,indices(i,:)+1]);
            end
            knots = knots(indices(i,:));
            equalKnotVecs = true;
            for jj = 1:d_p
                if any(nurbs1.number ~= nurbs2.number(indices(i,:))) || ...
                   any(nurbs1.degree ~= nurbs2.degree(indices(i,:))) || ...
                   any(abs(nurbs1.knots{jj} - knots{jj}) > Eps)
                    equalKnotVecs = false;
                    break
                end
            end
            isEqual = all(abs(nurbs1.coeffs(:) - temp(:)) < Eps) && equalKnotVecs;
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