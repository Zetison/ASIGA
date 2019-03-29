function nurbs = translateNURBS(nurbs,x0)

if isrow(x0)
    x0 = x0.';
end

nurbs.coeffs(1:3,:,:) = nurbs.coeffs(1:3,:,:) + repmat(x0,1,nurbs.number(1),nurbs.number(2));