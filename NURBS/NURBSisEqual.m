function isEqual = NURBSisEqual(nurbs1,nurbs2,purelyGeometric)
if nargin < 3
    purelyGeometric = false;
end
isEqual = nurbs1.d_p == nurbs2.d_p && ...
          nurbs1.d == nurbs1.d && ...
          all(sort(nurbs1.number) == sort(nurbs2.number)) && ...
          all(sort(nurbs1.degree) == sort(nurbs2.degree));
if ~isEqual
    return
end
Eps = 1e-10;
d_p = nurbs1.d_p;
if ~purelyGeometric
    for i = 1:d_p
        if any(abs(nurbs1.knots{i} - nurbs2.knots{i}) > Eps)
            isEqual = false;
        end
    end
    if ~isEqual
        return
    end
end
if purelyGeometric
    flips = cell(2^d_p,1);
    flips{1} = [];
    flips{2} = 1;
    flips{3} = 2;
    flips{4} = [1,2];
    flips{5} = 3;
    flips{6} = [1,3];
    flips{7} = [2,3];
    flips{8} = [1,2,3];
    indices = perms(1:d_p)+1;
    for i = 1:size(indices,1)
        for j = 1:2^d_p
            coeffs2 = nurbs2.coeffs;
            for jj = 1:numel(flips{j})
                coeffs2 = flip(coeffs2,flips{j}(jj)+1);
            end
            temp = permute(coeffs2,[1,indices(i,:)]);
            isEqual = all(abs(nurbs1.coeffs(:) - temp(:)) < Eps);
            if isEqual
                return
            end
        end
    end
else
    temp = abs(nurbs1.coeffs - nurbs2.coeffs) < Eps;
    isEqual = all(temp(:));
end