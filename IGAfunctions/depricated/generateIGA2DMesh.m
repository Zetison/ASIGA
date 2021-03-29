function [element, index, noElems] = generateIGA2DMesh(Xi, Eta, p, q, n, m)
% warning('This function is depricated. Use generateIGAmesh instead')

uniqueXiKnots   = unique(Xi);
uniqueEtaKnots   = unique(Eta);

noElemsXi       = length(uniqueXiKnots)-1; % # of elements xi dir.
noElemsEta       = length(uniqueEtaKnots)-1; % # of elements eta dir.

chan  = zeros(m,n);

count = 1;
for j = 1:m
    for i = 1:n
        chan(i,j) = count;
        count     = count + 1;
    end
end


% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeXi,elConnXi] = buildConnectivity(p,Xi,noElemsXi);
[elRangeEta,elConnEta] = buildConnectivity(q,Eta,noElemsEta);

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

index = zeros(noElems,2);

count = 1;
for j = 1:size(elRangeEta,1)
    for i = 1:size(elRangeXi,1)
        index(count,1) = i;
        index(count,2) = j;
        count = count + 1;
    end
end









