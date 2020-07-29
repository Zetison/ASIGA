function varCol = generateIGAmesh(varCol)


patches = varCol.patches;
switch patches{1}.nurbs.d_p
    case 1
        for patch = 1:numel(patches)
            Xi = patches{patch}.nurbs.knots{1};

            p = patches{patch}.nurbs.degree(1);   

            n = patches{patch}.nurbs.number(1);  
            
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

            patches{patch}.element = element;
            patches{patch}.index = index;
            patches{patch}.noElems = noElems;
            patches{patch}.elRange{1} = elRangeXi;
        end
    case 2
        for patch = 1:numel(patches)
            Xi = patches{patch}.nurbs.knots{1};
            Eta = patches{patch}.nurbs.knots{2};  

            p = patches{patch}.nurbs.degree(1);  
            q = patches{patch}.nurbs.degree(2);  

            n = patches{patch}.nurbs.number(1);  
            m = patches{patch}.nurbs.number(2);  

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
            counter = 1;
            for j = 1:m
                for i = 1:n
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

            index = zeros(noElems,2);
            counter = 1;

            for j = 1:noElemsEta
                for i = 1:noElemsXi
                    index(counter,1) = i;
                    index(counter,2) = j;
                    counter = counter + 1;
                end
            end

            patches{patch}.element = element;
            patches{patch}.index = index;
            patches{patch}.noElems = noElems;
            patches{patch}.elRange{1} = elRangeXi;
            patches{patch}.elRange{2} = elRangeEta;
        end
    case 3
        for patch = 1:numel(patches)
            Xi = patches{patch}.nurbs.knots{1};
            Eta = patches{patch}.nurbs.knots{2};
            Zeta = patches{patch}.nurbs.knots{3};   

            p = patches{patch}.nurbs.degree(1);  
            q = patches{patch}.nurbs.degree(2);  
            r = patches{patch}.nurbs.degree(3);  

            n = patches{patch}.nurbs.number(1);  
            m = patches{patch}.nurbs.number(2);  
            l = patches{patch}.nurbs.number(3);  

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

            counter = 1;
            for k = 1:l
                for j = 1:m
                    for i = 1:n
                        chan(i,j,k) = counter;
                        counter       = counter + 1;
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

            patches{patch}.element = element;
            patches{patch}.index = index;
            patches{patch}.noElems = noElems;
            patches{patch}.elRange{1} = elRangeXi;
            patches{patch}.elRange{2} = elRangeEta;
            patches{patch}.elRange{3} = elRangeZeta;
        end
end
varCol.patches = patches;





