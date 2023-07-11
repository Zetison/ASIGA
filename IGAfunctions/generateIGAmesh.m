function varCol = generateIGAmesh(varCol)


patches = varCol.patches;
for patch = 1:numel(patches)
    d_p = patches{patch}.nurbs.d_p;
    knots = patches{patch}.nurbs.knots;
    degree = patches{patch}.nurbs.degree;   
    number = patches{patch}.nurbs.number;  

    elRange = cell(1,d_p);
    elConn = cell(1,d_p);
    noElemsArr = zeros(1,d_p);
    for i = 1:d_p
        noElemsArr(i) = length(unique(knots{i}))-1;
        [elRange{i},elConn{i}] = buildConnectivity(degree(i),knots{i},noElemsArr(i));
    end

    noElems = prod(noElemsArr);
    element = zeros(noElems,prod(degree+1));
    chan = reshape(1:prod(number),[number,1]);
    index = zeros(noElems,d_p);
    
    switch d_p
        case 1
            e = 1;
            for ii = 1:noElemsArr(1)
                c = 1;
                uConn = elConn{1}(ii,:);

                for i = 1:length(uConn)
                    element(e,c) = chan(uConn(i));
                    c = c + 1;
                end
                e = e + 1;
            end

            index = (1:noElemsArr(1)).';
        case 2
            e = 1;
            for jj = 1:noElemsArr(2)
                vConn = elConn{2}(jj,:);
                for ii = 1:noElemsArr(1)
                    c = 1;
                    uConn = elConn{1}(ii,:);

                    for j = 1:length(vConn)
                        for i = 1:length(uConn)
                            element(e,c) = chan(uConn(i), vConn(j));
                            c = c + 1;
                        end
                    end
                    e = e + 1;
                end
            end

            [arr1,arr2] = ndgrid(1:noElemsArr(1),1:noElemsArr(2));
            index = [arr1(:),arr2(:)];
        case 3
            e = 1;
            for kk = 1:noElemsArr(3)
                wConn = elConn{3}(kk,:);
                for jj = 1:noElemsArr(2)
                    vConn = elConn{2}(jj,:);
                    for ii = 1:noElemsArr(1)
                        c = 1;
                        uConn = elConn{1}(ii,:);

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

            [arr1,arr2,arr3] = ndgrid(1:noElemsArr(1),1:noElemsArr(2),1:noElemsArr(3));
            index = [arr1(:),arr2(:),arr3(:)];
    end
    patches{patch}.element = element;
    patches{patch}.index = index;
    patches{patch}.noElems = noElems;
    patches{patch}.elRange = elRange;
end
varCol.patches = patches;





