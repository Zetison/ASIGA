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
            case 'xy'
                temp = coeffs(2,:,:,:);
                coeffs(2,:,:,:) = coeffs(1,:,:,:);
                coeffs(1,:,:,:) = temp;
        end
        nurbs{patch}.coeffs = coeffs;
    end
else
    coeffs = nurbs.coeffs;
    switch axis
        case 'x'
            coeffs(1,:,:) = -coeffs(1,:,:);
        case 'y'
            coeffs(2,:,:) = -coeffs(2,:,:);
        case 'z'
            coeffs(3,:,:) = -coeffs(3,:,:);
        case 'xy'
            temp = coeffs(2,:,:,:);
            coeffs(2,:,:,:) = coeffs(1,:,:,:);
            coeffs(1,:,:,:) = temp;
    end
    nurbs.coeffs = coeffs;
end