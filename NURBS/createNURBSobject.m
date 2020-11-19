function nurbs = createNURBSobject(coeffs,knots)
if nargin < 2
    knots = {};
end
if ~iscell(knots)
    knots = {knots};
end
np = size(coeffs);
d_p = length(knots);
number = np(2:end);

sizeKnots = zeros(1,d_p);
for i = 1:d_p % normalize knots
    knots{i} = (knots{i}-knots{i}(1))/(knots{i}(end)-knots{i}(1)); 
    sizeKnots(i) = numel(knots{i});
end
nurbs.d_p = d_p;
nurbs.d = np(1)-1;
nurbs.number = number; % array of the number of basis functions in each parametric direction
nurbs.degree = sizeKnots-number-1;
nurbs.knots = knots;
nurbs.coeffs = coeffs;
nurbs = {nurbs};