function newNurbs = elevateNURBSdegree(nurbs,degreeElevations)

d_p = numel(nurbs.knots);
d = size(nurbs.coeffs,1)-1;


coeffs = nurbs.coeffs;
coeffs(1:d,:,:,:) = coeffs(1:d,:,:,:).*repmat(coeffs(end,:,:,:),d,1,1,1);

knots = cell(1,d_p);
for i = 1:d_p
    if degreeElevations(i) == 0
        knots{i} = nurbs.knots{i};
    else   
        dimensions = size(coeffs);
        indices = 1:d_p+1;
        indices([i+1,d_p+1]) = [d_p+1,i+1];
        coeffs = permute(coeffs,indices);
        coeffs = reshape(coeffs,prod(dimensions(indices(1:end-1))),dimensions(i+1));
        [coeffs, knots{i}] = elevateBsplinesDegree(nurbs.number(i), nurbs.degree(i), nurbs.knots{i}, coeffs, degreeElevations(i));
        n = size(coeffs,2);
        coeffs = reshape(coeffs,[dimensions(indices(1:end-1)),n]);
        coeffs = permute(coeffs,indices);
    end
end
coeffs(1:d,:,:,:) = coeffs(1:d,:,:,:)./repmat(coeffs(end,:,:,:),d,1,1,1);

newNurbs = createNURBSobject(coeffs,knots); % construct new NURBS