function nurbsPatches = elevateNURBSdegree(nurbsPatches,degreeElevations)

for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    d = nurbs.d;
    coeffs = nurbs.coeffs;
    d_p = nurbs.d_p;
    d_pp1 = d_p+1;

    coeffs = NURBSprojection(coeffs,d,d_p,'unproject');

    knots = cell(1,d_p);
    if numel(degreeElevations) < d_p
        degreeElevs = zeros(1,d_p);
        tmp = numel(degreeElevations);
        degreeElevs(1:tmp) = degreeElevations;
        degreeElevs(tmp+1:d_p) = degreeElevations(end);
    else
        degreeElevs = degreeElevations;
    end
    for i = 1:d_p
        if degreeElevs(i) == 0
            knots{i} = nurbs.knots{i};
        else   
            dimensions = size(coeffs);
            indices = 1:d_pp1;
            indices([i+1,d_pp1]) = [d_pp1,i+1];
            coeffs = permute(coeffs,indices);
            coeffs = reshape(coeffs,prod(dimensions(indices(1:end-1))),dimensions(i+1));
            [coeffs, knots{i}] = elevateBsplinesDegree(nurbs.number(i), nurbs.degree(i), nurbs.knots{i}, coeffs, degreeElevs(i));
            dimensions(i+1) = size(coeffs,2);
            coeffs = reshape(coeffs,dimensions(indices));
            coeffs = permute(coeffs,indices);
        end
    end
    nurbs.coeffs = coeffs;
    coeffs = NURBSprojection(coeffs,d,d_p,'project');

    nurbsPatches{patch}.coeffs = coeffs;
    nurbsPatches{patch}.knots = knots;
    np = size(coeffs);
    nurbsPatches{patch}.number = np(2:end);
    nurbsPatches{patch}.degree = nurbsPatches{patch}.degree + degreeElevs;
end


