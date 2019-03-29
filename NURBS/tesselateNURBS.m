function [X, elements,nurbs] = tesselateNURBS(nurbs)



nurbs = repeatKnots(nurbs,'linear_FEM');
nurbs = degenerateIGAtoFEM_Grev(nurbs,'linear_FEM');

Xi = nurbs.knots{1};
Eta = nurbs.knots{2};

noXiKnots = length(unique(Xi));
noEtaKnots = length(unique(Eta));

noVisElems  = (noXiKnots-1)*(noEtaKnots-1);
elements = zeros(noVisElems,4);
eVis = 1;

for j = 1:noEtaKnots-1
    for i = 1:noXiKnots-1
        elements(eVis,1) = i   +   (j-1)*noXiKnots;
        elements(eVis,2) = i+1 +   (j-1)*noXiKnots;
        elements(eVis,3) = i+1 +       j*noXiKnots;
        elements(eVis,4) = i   +       j*noXiKnots;

        eVis = eVis + 1;
    end
end

p = nurbs.degree(1);
q = nurbs.degree(2);

n = length(Xi)-p-1;
m = length(Eta)-q-1;
X = zeros(n*m,3);

count = 0;
for j=1:m
    X(n*count+1:n*(count+1),:) = nurbs.coeffs(1:3,:,j)';
    count = count+1;
end