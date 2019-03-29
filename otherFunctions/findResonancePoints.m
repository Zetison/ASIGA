function varCol = findResonancePoints(varCol, sp, d, elementList)
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
index = varCol.index;
nurbs = varCol.nurbs;
cps.cp_a = zeros(2,2);
cps.cp_b = zeros(2,2);
cps.cp_c = zeros(2,2);
cps.cp_d = zeros(2,2);
if nargin < 4
    elementList = 1:noElems;
end
for e = elementList
    idXi = index(e,1);
    idEta = index(e,2);
    
    xiE = elRangeXi(idXi,:);
    etaE = elRangeEta(idEta,:);
    
    
    maxItrs = 100;
    Eps = 1e4*eps;
    for i = 1:2
        x = xiE(i);
        dgdy = @(y) sphereMappingNURBS(x, y, nurbs, sp, d, 'dgdy');
        dgdy2 = @(y) sphereMappingNURBS(x, y, nurbs, sp, d, 'dgdy2');
        
        y_res = newtonsMethod(dgdy,dgdy2,sum(etaE)/2,maxItrs,Eps);
        if dgdy(y_res) > Eps
            y_res = Inf;
        end
        
        if abs(dgdy2(y_res)) < eps
            r = 2;
        else
            r = 1;
        end
        if i == 1
            cps(e).cp_a = [etaE(1), etaE(2);
                            0, 0];
            if etaE(1) < y_res && y_res < etaE(2)
                cps(e).cp_a = [etaE(1), y_res, etaE(2);
                               0, r, 0];
            elseif y_res == etaE(2)
                cps(e).cp_a(2,:) = [etaE(2); r];
            elseif y_res == etaE(1)
                cps(e).cp_a(1,:) = [etaE(1); r];
            end
        else
            cps(e).cp_b = [etaE(1), etaE(2);
                            0, 0];
            if etaE(1) < y_res && y_res < etaE(2)
                cps(e).cp_b = [etaE(1), y_res, etaE(2);
                               0, r, 0];
            elseif y_res == etaE(2)
                cps(e).cp_b(2,:) = [etaE(2); r];
            elseif y_res == etaE(1)
                cps(e).cp_b(1,:) = [etaE(1); r];
            end
        end
    end
    for i = 1:2
        y = etaE(i);
        dgdx = @(x) sphereMappingNURBS(x, y, nurbs, sp, d, 'dgdx');
        dgdx2 = @(x) sphereMappingNURBS(x, y, nurbs, sp, d, 'dgdx2');
        
        x_res = newtonsMethod(dgdx,dgdx2,sum(xiE)/2,maxItrs,Eps);
        if dgdx(x_res) > Eps
            x_res = Inf;
        end
        
        if abs(dgdx2(x_res)) < eps
            r = 2;
        else
            r = 1;
        end
        if i == 1
            cps(e).cp_c = [xiE(1), xiE(2);
                            0, 0];
            if xiE(1) < x_res && x_res < xiE(2)
                cps(e).cp_c = [xiE(1), x_res, xiE(2);
                               0, r, 0];
            elseif x_res == xiE(1)
                cps(e).cp_c(2,:) = [xiE(2); r];
            elseif x_res == xiE(2)
                cps(e).cp_c(1,:) = [xiE(1); r];
            end
        else
            cps(e).cp_d = [xiE(1), xiE(2);
                            0, 0];
            if xiE(1) < x_res && x_res < xiE(2)
                cps(e).cp_d = [xiE(1), x_res, xiE(2);
                               0, r, 0];
            elseif x_res == xiE(1)
                cps(e).cp_d(2,:) = [xiE(2); r];
            elseif x_res == xiE(2)
                cps(e).cp_d(1,:) = [xiE(1); r];
            end
        end
    end
end
% if numel(elementList) == 1
%     cps = cps(e);
% end
varCol.cps = cps;