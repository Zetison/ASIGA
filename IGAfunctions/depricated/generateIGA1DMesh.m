function varCol = generateIGA1DMesh(varCol)
error('Depricated. Use generateIGAmesh instead')
Xi = varCol.patches{1}.nurbs.knots;
p = varCol.patches{1}.nurbs.degree;
n = varCol.patches{1}.nurbs.number;
uniqueXiKnots   = unique(Xi);

noElems       = length(uniqueXiKnots)-1; % # of elements xi dir.

chan = 1:n;

% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeXi,elConnXi] = buildConnectivity(p,Xi,noElems);

element = zeros(noElems,p+1);

e = 1;
for ii = 1:noElems
    c = 1;
    uConn = elConnXi(ii,:);

    for i = 1:length(uConn)
        element(e,c) = chan(uConn(i));
        c = c + 1;
    end
    e = e + 1;
end

index = zeros(noElems,2);

count = 1;
for i = 1:size(elRangeXi,1)
    index(count,1) = i;
    count = count + 1;
end

varCol.element = element;
varCol.index = index;
varCol.noElems = noElems;
varCol.elRangeXi = elRangeXi;





