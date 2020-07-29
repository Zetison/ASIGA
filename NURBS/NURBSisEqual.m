function isEqual = NURBSisEqual(nurbs1,nurbs2,purelyGeometric)
if nargin < 3
    purelyGeometric = false;
end
isEqual = nurbs1.d_p == nurbs2.d_p && ...
          nurbs1.d == nurbs1.d && ...
          all(nurbs1.degree == nurbs2.degree) && ...
          all(nurbs1.number == nurbs2.number);
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
    indices = perms(1:d_p)+1;
    for i = 1:size(indices,1)
        for j = 2:d_p+2
            temp = abs(nurbs1.coeffs - permute(flip(nurbs2.coeffs,j),[1,indices(i,:)])) < Eps;
            isEqual = all(temp(:));
            if isEqual
                return
            end
        end
    end
else
    temp = abs(nurbs1.coeffs - nurbs2.coeffs) < Eps;
    isEqual = all(temp(:));
end