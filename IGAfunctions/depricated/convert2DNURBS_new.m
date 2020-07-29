function varCol = convert2DNURBS_new(nurbs, varCol)
error('Depricated, use convertNURBS instead')

p = nurbs.degree(1);
q = nurbs.degree(2);

Xi   = nurbs.knots{1};
Eta  = nurbs.knots{2};

n = length(Xi)-p-1;
m = length(Eta)-q-1;

weights = reshape(nurbs.coeffs(3,:,:),n*m,1);

controlPts = zeros(n*m,2);

count = 0;
for j=1:m
    controlPts(n*count+1:n*(count+1),:) = nurbs.coeffs(1:2,:,j)';
    count = count+1;
end

varCol.weights = weights;
varCol.controlPts = controlPts;


varCol.nurbs = nurbs;
varCol.noCtrlPts = nurbs.number(1)*nurbs.number(2);
varCol.noDofs = varCol.dimension*varCol.noCtrlPts;