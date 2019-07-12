function nurbs = translateNURBS(nurbs,x0)

if isrow(x0)
    x0 = x0.';
end

if ~iscell(nurbs)
    nurbs = {nurbs};
end
for patch = 1:numel(nurbs)
    coeffs = nurbs{patch}.coeffs;
    number = nurbs{patch}.number;
    switch nurbs{patch}.type
        case '3Dsurface'
            coeffs(1:3,:,:) = coeffs(1:3,:,:) + repmat(x0,1,number(1),number(2));
        case '3Dvolume'
            coeffs(1:3,:,:,:) = coeffs(1:3,:,:,:) + repmat(x0,1,number(1),number(2),number(3));
    end
    nurbs{patch}.coeffs = coeffs;
end