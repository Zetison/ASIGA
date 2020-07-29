function nurbsPatches = elevateNURBSdegree(nurbsPatches,degreeElevations)

for patch = 1:numel(nurbsPatches)
    nurbs = nurbsPatches{patch};
    d = nurbs.d;
    d_p = nurbs.d_p;
    d_pp1 = d_p+1;

    coeffs = nurbs.coeffs;
    coeffs = subasgnArr(coeffs,slc(coeffs,1:d).*repmat(slc(coeffs,d+1),[d,ones(1,d_p)]),1:d);

    knots = cell(1,d_p);
    for i = 1:d_p
        if degreeElevations(i) == 0
            knots{i} = nurbs.knots{i};
        else   
            dimensions = size(coeffs);
            indices = 1:d_pp1;
            indices([i+1,d_pp1]) = [d_pp1,i+1];
            coeffs = permute(coeffs,indices);
            coeffs = reshape(coeffs,prod(dimensions(indices(1:end-1))),dimensions(i+1));
            [coeffs, knots{i}] = elevateBsplinesDegree(nurbs.number(i), nurbs.degree(i), nurbs.knots{i}, coeffs, degreeElevations(i));
            dimensions(i+1) = size(coeffs,2);
            coeffs = reshape(coeffs,dimensions(indices));
            coeffs = permute(coeffs,indices);
        end
    end
    coeffs = subasgnArr(coeffs,slc(coeffs,1:d)./repmat(slc(coeffs,d+1),[d,ones(1,d_p)]),1:d);

    nurbsPatches(patch) = createNURBSobject(coeffs,knots);
end


