function isTrue = NURBSisSingular(nurbs)
error('Depricated, use NURBSisDegenerate() instead')
d = nurbs.d;
coeffs = slc(nurbs.coeffs,1:d,1);
sizes = size(coeffs);
isTrue = false;
for i = 1:d_p
    dims = ones(1,d_p);
    dims(i) = sizes(i);
    if all(repmat(slc(coeffs,1,2),[1,dims])-coeffs < 1e-10)
        isTrue = true;
        return
    end
end