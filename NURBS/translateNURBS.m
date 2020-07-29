function nurbs = translateNURBS(nurbs,x0)

if isrow(x0)
    x0 = x0.';
end
if numel(x0) == 3
    nurbs = ensure3DNURBS(nurbs);
end

if ~iscell(nurbs)
    nurbs = {nurbs};
end
for patch = 1:numel(nurbs)
    coeffs = nurbs{patch}.coeffs;
    number = nurbs{patch}.number;
    d = size(coeffs,1)-1;
    
    coeffs = subasgnArr(coeffs,slc(coeffs, 1:d) + repmat(x0,[1,number]),1:d);
    nurbs{patch}.coeffs = coeffs;
end