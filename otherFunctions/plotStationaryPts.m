function plotStationaryPts(varCol, sp, d, e, x_arr)

elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
index = varCol.index;
nurbs = varCol.nurbs;

idXi = index(e,1);
idEta = index(e,2);

xiE = elRangeXi(idXi,:);
etaE = elRangeEta(idEta,:);


maxItrs = 100;
Eps = 1e4*eps;
for j = 1:numel(x_arr)
    x = x_arr(j);
    dgdy = @(y) sphereMappingNURBS(x, y, nurbs, sp, d, xiE, etaE, 'dgdy');
    dgdy2 = @(y) sphereMappingNURBS(x, y, nurbs, sp, d, xiE, etaE, 'dgdy2');

    y_res = newtonsMethod(dgdy,dgdy2,0,maxItrs,Eps);
    
    if -1 < y_res && y_res < 1
        plot(x,y_res,'*')
    end
end
y_arr = x_arr;
for j = 1:numel(y_arr)
    y = y_arr(j);
    dgdx = @(x) sphereMappingNURBS(x, y, nurbs, sp, d, xiE, etaE, 'dgdx');
    dgdx2 = @(x) sphereMappingNURBS(x, y, nurbs, sp, d, xiE, etaE, 'dgdx2');

    x_res = newtonsMethod(dgdx,dgdx2,0,maxItrs,Eps);
    
    if -1 < x_res && x_res < 1
        plot(x_res,y,'*')
    end
end