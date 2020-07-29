function [element, index, noElems, ...
          elRangeXi, elRangeEta, elRangeZeta, childrenNodes, ...
          elConnXi, elConnEta, elConnZeta, ...
          noElemsXi, noElemsEta, noElemsZeta] = generateIGA3DMeshSpecial(Xi, Eta, Zeta, p, q, r, l, m, n)
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
childrenNodes = zeros(1,2*(n-1)*l+(m-2)*l);

counter = 1;
counter2 = 1;
for k = 1:l
    for j = 1:m
        for i = 1:n
            if j == 1 || j == m
                chan(:,j,k) = counter;
                childrenNodes(counter2:counter2+(n-2)) = (counter+1):(counter+n-1);
                counter = counter + n;
                counter2 = counter2 + n-1;
                break
            elseif i == 1
                chan([1 n],j,k) = counter;
            elseif i == n
                childrenNodes(counter2) = counter;
                counter2 = counter2 + 1;
                counter = counter + 1;
                continue
            end
            
            chan(i,j,k) = counter;
            counter     = counter + 1;
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
counter = 1;

for k = 1:noElemsZeta
    for j = 1:noElemsEta
        for i = 1:noElemsXi
            index(counter,1) = i;
            index(counter,2) = j;
            index(counter,3) = k;
            counter = counter + 1;
        end
    end
end






