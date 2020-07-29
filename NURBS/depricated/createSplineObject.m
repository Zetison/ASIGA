function splineObj = createSplineObject(coeffs,knots)
error('Depricated, use createNURBSobject instead')
np = size(coeffs);
if iscell(knots)
    % Construct spline volume in 3D
    splineObj.type = '3Dvolume';
    n = np(2);
    m = np(3);
    l = np(4);
    splineObj.coeffs = coeffs;
    Xi = knots{1};
    Eta = knots{2};
    Zeta = knots{3};    
    splineObj.knots = {Xi/Xi(end) Eta/Eta(end) Zeta/Zeta(end)};
    p = size(knots{1},2)-n-1;
    q = size(knots{2},2)-m-1;
    r = size(knots{3},2)-l-1;
    splineObj.degree = [p q r];
    splineObj.number = [n m l];
else
    % constructing a spline curve in 2D
    splineObj.type = '2Dcurve';
    splineObj.coeffs = coeffs;
    n = np(2);
    splineObj.number = n;
    splineObj.degree = size(knots,2)-n-1;
    splineObj.knots = knots/knots(end);
end
