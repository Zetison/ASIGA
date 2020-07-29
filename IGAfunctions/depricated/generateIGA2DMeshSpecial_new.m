function varCol = generateIGA2DMeshSpecial_new(varCol)
error('Depricated. Use generateIGAmesh instead')

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};

p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);

n = varCol.nurbs.number(1);  
m = varCol.nurbs.number(2); 

uniqueXiKnots   = unique(Xi);
uniqueEtaKnots   = unique(Eta);

noElemsXi       = length(uniqueXiKnots)-1; % # of elements xi dir.
noElemsEta       = length(uniqueEtaKnots)-1; % # of elements eta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points


chan  = zeros(n,m);
childrenNodes = zeros(1,2*(n-1)+m-2);

counter = 1;
counter2 = 1;
for j = 1:m
    for i = 1:n
        if j == 1 || j == m
            chan(:,j) = counter;
            childrenNodes(counter2:counter2+(n-2)) = (counter+1):(counter+n-1);
            counter = counter + n;
            counter2 = counter2 + n-1;
            break
        elseif i == 1
            chan([1 n],j) = counter;
        elseif i == n
            childrenNodes(counter2) = counter;
            counter2 = counter2 + 1;
            counter = counter + 1;
            continue
        end

        chan(i,j) = counter;
        counter     = counter + 1;
    end
end



% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeXi,elConnXi] = buildConnectivity(p,Xi,noElemsXi);
[elRangeEta,elConnEta] = buildConnectivity(q,Eta,noElemsEta);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

noElems = noElemsXi * noElemsEta;
element = zeros(noElems,(p+1)*(q+1));

e = 1;
for jj = 1:noElemsEta
    vConn = elConnEta(jj,:);
    for ii = 1:noElemsXi
        c = 1;
        uConn = elConnXi(ii,:);

        for j = 1:length(vConn)
            for i = 1:length(uConn)
                element(e,c) = chan(uConn(i), vConn(j));
                c = c + 1;
            end
        end
        e = e + 1;
    end
end

index = zeros(noElems,3);
counter = 1;

for j = 1:noElemsEta
    for i = 1:noElemsXi
        index(counter,1) = i;
        index(counter,2) = j;
        counter = counter + 1;
    end
end

varCol.element = element;
varCol.index = index;
varCol.noElems = noElems;
varCol.elRangeXi = elRangeXi;
varCol.elRangeEta = elRangeEta;
varCol.childrenNodes = childrenNodes;





