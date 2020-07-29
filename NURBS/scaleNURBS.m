function nurbs = scaleNURBS(nurbs,c)

if all(c == 1)
    return
end
if isrow(c)
    c = c.';
end
for patch = 1:numel(nurbs)
    coeffs = nurbs{patch}.coeffs;
    d = nurbs{patch}.d;
    if numel(c) == 1
        cVec = c*ones(d,1);
    else
        cVec = c;
    end
    coeffs = subasgnArr(coeffs,cVec.*slc(coeffs, 1:d),1:d);

    nurbs{patch}.coeffs = coeffs;
end