function varargout = evaluateNURBS(nurbs,xi,n,U)
% This routine evaluates all derivatives up to order n for a single NURBS patch at the parametric points xi (having sizes (npts, d_p))
% If U is provided, it is used as coefficients instead of the NURBS coefficients.

if nargin < 3
    n = 0;
end
d = nurbs.d;
d_p = nurbs.d_p;
knots = nurbs.knots;
degree = nurbs.degree;
number = nurbs.number;
I = zeros(size(xi,1),d_p);
for i = 1:d_p
    I(:,i) = findKnotSpan(number(i), degree(i), xi(:,i), knots{i});
end
weights = reshape(nurbs.coeffs(d+1,:),1,[]);

d_p = numel(degree);
order = degree+1;
n_en = prod(order);
npts = size(xi,1);

w_i = ones([npts,order]);
shift = 1;
for i = 1:d_p
    temp_i = ones(1,d_p);
    temp_i(i) = degree(i)+1;
    w_i = w_i + (reshape(I(:,i)-degree(i)+(0:degree(i)),[npts,temp_i])-1)*shift;
    shift = shift*number(i);
end
w = weights(w_i);

RdR = NURBSbasis(I,xi,degree,knots,w,n);
if nargin < 4
    coeffs = nurbs.coeffs(1:d,:).';
    U = reshape(coeffs(w_i,:),npts,n_en,d);
else
    U = reshape(U(w_i,:),npts,n_en,size(U,2));
end
d_p = numel(knots);
varargout = cell(1,d_p+1);
varargout{1} = reshape(sum(RdR{1}.*U,2),[],size(U,3));
for i = 1:d_p
    for j = 1:n
        varargout{i+1}(:,:,j) = reshape(sum(RdR{i+1}(:,:,j).*U,2),[],size(U,3));
    end
end
