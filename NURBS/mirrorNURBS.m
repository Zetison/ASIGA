function nurbs = mirrorNURBS(nurbs,axis)

switch axis
    case 'x'
        nurbs.coeffs(1,:,:) = -nurbs.coeffs(1,:,:);
    case 'y'
        nurbs.coeffs(2,:,:) = -nurbs.coeffs(2,:,:);
    case 'z'
        nurbs.coeffs(3,:,:) = -nurbs.coeffs(3,:,:);
end