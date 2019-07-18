function nurbs = mirrorNURBS(nurbs,axis)

if iscell(nurbs)
    for patch = 1:numel(nurbs)
        coeffs = nurbs{patch}.coeffs;
        switch axis
            case 'x'
                coeffs(1,:,:,:) = -coeffs(1,:,:,:);
            case 'y'
                coeffs(2,:,:,:) = -coeffs(2,:,:,:);
            case 'z'
                coeffs(3,:,:,:) = -coeffs(3,:,:,:);
        end
        nurbs{patch}.coeffs = coeffs;
    end
else
    switch axis
        case 'x'
            nurbs.coeffs(1,:,:) = -nurbs.coeffs(1,:,:);
        case 'y'
            nurbs.coeffs(2,:,:) = -nurbs.coeffs(2,:,:);
        case 'z'
            nurbs.coeffs(3,:,:) = -nurbs.coeffs(3,:,:);
    end
end