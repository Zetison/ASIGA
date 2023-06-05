function nurbsPatches = insertKnotsInNURBS(nurbsPatches,newKnots)
noPatches = numel(nurbsPatches);
uniformRef = false;
if iscell(newKnots) && ~iscell(newKnots{1})
    temp = newKnots;
    newKnots = cell(1,noPatches);
    [newKnots{1:noPatches}] = deal(temp);
elseif ~iscell(newKnots)
    uniformRef = true;
end

for patch = 1:noPatches
    nurbs = nurbsPatches{patch};
    d = nurbs.d;
    coeffs = nurbs.coeffs;
    number = nurbs.number;
    d_p = nurbs.d_p;
    d_pp1 = d_p+1;
    coeffs = NURBSprojection(coeffs,d,d_p,'unproject');
    knots = cell(1,d_p);
    for i = 1:d_p
        if uniformRef
            if numel(newKnots) < d_p
                newKnots(numel(newKnots):d_p) = newKnots(end);
            end
            newKnotsPatch = insertUniform2(nurbs.knots{i}, newKnots(i));
        else
            if isempty(newKnots{patch})
                newKnotsPatch = [];
            else
                newKnotsPatch = newKnots{patch}{i};
            end
        end
        if isempty(newKnotsPatch)
            knots{i} = nurbs.knots{i};
        else   
            dimensions = size(coeffs);
            indices = 1:d_pp1;
            indices([i+1,d_pp1]) = [d_pp1,i+1];
            coeffs = permute(coeffs,indices);
            coeffs = reshape(coeffs,prod(dimensions(indices(1:end-1))),number(i));
            [coeffs,knots{i}] = insertKnotsInBsplines(number(i), nurbs.degree(i), nurbs.knots{i}, newKnotsPatch, coeffs);
            dimensions(i+1) = size(coeffs,2);
            coeffs = reshape(coeffs,dimensions(indices));
            coeffs = permute(coeffs,indices);
        end
    end
    coeffs = NURBSprojection(coeffs,d,d_p,'project');
    
    nurbsPatches{patch}.coeffs = coeffs;
    nurbsPatches{patch}.knots = knots;
    np = size(coeffs);
    nurbsPatches{patch}.number = np(2:end);
end










