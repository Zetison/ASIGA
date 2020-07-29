function varargout = evaluateNURBS(nurbs, parm_pt,n)

if nargin < 3
    n = 0;
end
d = nurbs.d;
d_p = nurbs.d_p;
knots = nurbs.knots;
degree = nurbs.degree;
number = nurbs.number;
coeffs = nurbs.coeffs;
I = zeros(1,d_p);
ix = cell(1,d_p+1);
for i = 1:d_p
    I(i) = findKnotSpan(number(i), degree(i), parm_pt(1,i), knots{i});
    ix{i+1} = I(i) - degree(i) + (1:degree(i)+1) - 1;
end
nen = prod(degree+1);
ix{1} = 1:d;
P = slc(coeffs,ix,1:d_p+1);
P = reshape(P,d,nen);
ix{1} = d+1;
weights = slc(coeffs,ix,1:d_p+1);
weights = reshape(weights,1,nen);
varargout = evalNURBSloc(I,parm_pt,degree,knots,weights.',P.',n);
