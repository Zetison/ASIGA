function gk_arr = findCGpoints(varCol, U)
p = varCol.nurbs.degree;
noElems = varCol.noElems;
k = varCol.k;

elRangeXi = varCol.elRangeXi;
index = varCol.index;

nppe = 2*p-1;
g0_arr = zeros(noElems,nppe);
for e = 1:noElems
    idXi = index(e);

    Xi_e = elRangeXi(idXi,:);

    g0_arr(e,:) = linspace2(Xi_e(1),Xi_e(2),nppe);
end  
g0_arr = reshape(g0_arr,noElems*nppe,1);

gk_arr = zeros(size(g0_arr));
parfor i = 1:length(g0_arr)
    gk = g0_arr(i);
    gk_prev = inf;
    counter = 1;
    while abs(gk_prev-gk) > 10*eps && counter < 200
        gk_prev = gk;
        if gk <= 0
            gk = eps;
        elseif gk > 1
            gk = 1-eps;
        end
        [u, v, dudx, d2udx2, d3udx3] = numericalSolEval1D(gk, varCol, U);
        R_h = d2udx2 + k^2*u + varCol.f(v);
        dR_hdx = d3udx3 + k^2*dudx + varCol.dfdx(v);
        
        gk = gk - R_h/dR_hdx;
        counter = counter + 1;
    end
    gk_arr(i) = gk;
end
gk_arr = sort(uniquetol(gk_arr,1e-3/noElems/nppe));