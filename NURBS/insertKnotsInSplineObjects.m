function newSplineObj = insertKnotsInSplineObjects(splineObj,newKnots)

if strcmp(splineObj.type, '3Dvolume')
    % NURBS represents a volume
    n = splineObj.number(1);
    m = splineObj.number(2);
    l = splineObj.number(3);

    coeffs = splineObj.coeffs;
    % Insert knot in xi direction
    if isempty(newKnots{1})
        knots{1} = splineObj.knots{1};
    else   
        coeffs = permute(coeffs,[1 3 4 2]);
        coeffs = reshape(coeffs,3*m*l,n);
        [coeffs,knots{1}] = insertKnotsInBsplines(splineObj.number(1), splineObj.order(1), splineObj.knots{1}, newKnots{1}, coeffs);
        n = size(coeffs,2);
        coeffs = reshape(coeffs,[3 m l n]);
        coeffs = permute(coeffs,[1 4 2 3]);
    end
    % Insert knot in eta direction
    if isempty(newKnots{2})
        knots{2} = splineObj.knots{2};
    else   
        coeffs = permute(coeffs,[1 2 4 3]);
        coeffs = reshape(coeffs,3*n*l,m);
        [coeffs,knots{2}] = insertKnotsInBsplines(splineObj.number(2), splineObj.order(2), splineObj.knots{2}, newKnots{2}, coeffs);
        m = size(coeffs,2);
        coeffs = reshape(coeffs,[3 n l m]);
        coeffs = permute(coeffs,[1 2 4 3]);
    end
    % Insert knot in zeta direction
    if isempty(newKnots{3})
        knots{3} = splineObj.knots{3};
    else   
        coeffs = reshape(coeffs,3*n*m,l);
        [coeffs,knots{3}] = insertKnotsInBsplines(splineObj.number(3), splineObj.order(3), splineObj.knots{3}, newKnots{3}, coeffs);
        l = size(coeffs,2);
        coeffs = reshape(coeffs,[3 n m l]);
    end
elseif strcmp(splineObj.type, '2Dcurve')
    % NURBS represents a curve in 2D
    if isempty(newKnots)
        coeffs = splineObj.coefs;
        knots = splineObj.knots;
    else
        coeffs = splineObj.coeffs';
        [coeffs,knots] = insertKnotsInBsplines(splineObj.number, splineObj.order, splineObj.knots, newKnots, coeffs);
        coeffs = coeffs';
    end
end

% construct new NURBS
newSplineObj = createSplineObject(coeffs,knots);










