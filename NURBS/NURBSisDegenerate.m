function isDegenerate = NURBSisDegenerate(nurbs)

d = nurbs.d;
d_p = nurbs.d_p;
coeffs = slc(nurbs.coeffs,1:d,1);
sizes = size(coeffs);
isDegenerate = false;
for i = 1:d_p
    dims = ones(1,d_p);
    dims(i) = sizes(i);
    diffs = abs(repmat(slc(coeffs,1,i+1),[1,dims]) - coeffs) < 1e-10;
    if all(diffs(:))
        isDegenerate = true;
        return
    end
end