function [element, index, noElems, ...
          elRangeXi, elRangeEta, elRangeZeta, ...
          elConnXi, elConnEta, elConnZeta, ...
          noElemsXi, noElemsEta, noElemsZeta] = generateIGA3DMesh(Xi, Eta, Zeta, p, q, r, l, m, n)
error('Depricated. Use generateIGAmesh instead')

uniqueXiKnots   = unique(Xi);
uniqueEtaKnots   = unique(Eta);
uniqueZetaKnots   = unique(Zeta);

noElemsXi       = length(uniqueXiKnots)-1; % # of elements xi dir.
noElemsEta       = length(uniqueEtaKnots)-1; % # of elements eta dir.
noElemsZeta       = length(uniqueZetaKnots)-1; % # of elements zeta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points


chan  = zeros(n,m,l);

count = 1;
for k = 1:l
    for j = 1:m
        for i = 1:n
            chan(i,j,k) = count;
            count       = count + 1;
        end
    end
end


% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeXi,elConnXi] = buildConnectivity(p,Xi,noElemsXi);
[elRangeEta,elConnEta] = buildConnectivity(q,Eta,noElemsEta);
[elRangeZeta,elConnZeta] = buildConnectivity(r,Zeta,noElemsZeta);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

noElems = noElemsXi * noElemsEta * noElemsZeta;
element = zeros(noElems,(p+1)*(q+1)*(r+1));

e = 1;
for kk = 1:noElemsZeta
    wConn = elConnZeta(kk,:);
    for jj = 1:noElemsEta
        vConn = elConnEta(jj,:);
        for ii = 1:noElemsXi
            c = 1;
            uConn = elConnXi(ii,:);
            
            for k = 1:length(wConn)
                for j = 1:length(vConn)
                    for i = 1:length(uConn)
                        element(e,c) = chan(uConn(i), vConn(j), wConn(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

index = zeros(noElems,3);
count = 1;

for k = 1:noElemsZeta
    for j = 1:noElemsEta
        for i = 1:noElemsXi
            index(count,1) = i;
            index(count,2) = j;
            index(count,3) = k;
            count = count + 1;
        end
    end
end






