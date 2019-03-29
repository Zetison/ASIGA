function nurbs = getBCInternalWaterData(solid,alignWithAxis)
if nargin < 2
    alignWithAxis = 'Xaxis';
end

Xi = solid.knots{1};
Eta = solid.knots{2};
Zeta = [0 0 1 1];

coeffs_solid = solid.coeffs;
coeffs = zeros(4,size(coeffs_solid,2),size(coeffs_solid,3),2);
coeffs(:,:,:,2) = coeffs_solid(:,:,:,1);

coeffs([1,4],:,:,1) = coeffs([1,4],:,:,2);
h = abs(coeffs_solid(3,1,3,1));

coeffs(1,:,1:3,1) = coeffs(1,:,1:3,1)+h;
coeffs(1,:,end-1,1) = coeffs(1,:,end-2,1);
coeffs(1,:,end,1) = coeffs(1,:,end-2,1);

switch alignWithAxis
    case 'Xaxis' 
        % Nothing to be done
    case 'Yaxis' 
        temp = coeffs(1,:,:,:);
        coeffs(1,:,:,:) = coeffs(2,:,:,:);
        coeffs(2,:,:,:) = coeffs(3,:,:,:);
        coeffs(3,:,:,:) = temp;
    case 'Zaxis'
        temp = coeffs(1,:,:,:);
        coeffs(1,:,:,:) = coeffs(2,:,:,:);
        coeffs(2,:,:,:) = coeffs(3,:,:,:);
        coeffs(3,:,:,:) = temp;
end

nurbs = createNURBSobject(coeffs,{Xi, Eta, Zeta});
