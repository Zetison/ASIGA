function vargout = evaluateNURBSvec(nurbs, parm_pt,n)

if nargin < 3
    n = 0;
end
d = nurbs.d;
d_p = nurbs.d_p;
knots = nurbs.knots;
degree = nurbs.degree;
number = nurbs.number;
coeffs = nurbs.coeffs;
ix = cell(1,d_p+1);
tic
I = zeros(size(parm_pt,1),d_p);
I2 = zeros(size(parm_pt,1),d_p);
for i = 1:d_p
    for j = 1:size(parm_pt,1)
        I(j,i) = findKnotSpan(number(i), degree(i), parm_pt(j,i), knots{i});
    end
    I2(:,i) = findKnotSpanVec(number(i), degree(i), parm_pt(:,i), knots{i});
    ix{i+1} = I(:,i) - degree(i) + (1:degree(i)+1) - 1;
end
n_en = prod(degree+1);
ix{1} = 1:d;
P = slc(coeffs,ix,1:d_p+1);
P = reshape(P,d,n_en);
ix{1} = d+1;
weights = slc(coeffs,ix,1:d_p+1);
weights = reshape(weights,1,n_en);

i = 1;
tic
B = zeros(size(parm_pt,1),3,2);
for j = 1:size(parm_pt,1)
    B(j,:,:) = BsplineBasis(I(j,i), parm_pt(j,i), degree(i), knots{i}, n);
end
toc
tic
B2 = BsplineBasisVec(I2(:,i), parm_pt(:,i), degree(i), knots{i}, n);
toc



R = NURBSbasis(I,xi,degree,knots,weights,n);
d_p = numel(knots);
vargout = cell(1,d_p+1);
vargout{1} = R{1}*U;
for i = 1:d_p
    for j = 1:n
        vargout{i+1}(:,:,j) = R{i+1}(:,:,j)*U;
    end
end
