function nurbs = rotateNURBS(nurbs,theta,rotAxis)

if ~iscell(nurbs)
    nurbs = {nurbs};
end
for patch = 1:numel(nurbs)
    coeffs = nurbs{patch}.coeffs;
    R_x = rotationMatrix(theta, rotAxis);
    switch nurbs{patch}.type
        case '3Dsurface'
            for j = 1:size(coeffs,3)
                coeffs(1:3,:,j) = R_x*coeffs(1:3,:,j);
            end
        case '3Dvolume'
            for i = 1:size(coeffs,3)
                for j = 1:size(coeffs,4)
                    coeffs(1:3,:,i,j) = R_x*coeffs(1:3,:,i,j);
                end
            end
    end

    nurbs{patch}.coeffs = coeffs;
end