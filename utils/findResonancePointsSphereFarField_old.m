function cps = findResonancePointsSphereFarField_old(dmx, R, etaShift, elementList, xiE, etaE)
cps.cp_a = zeros(2,2);
cps.cp_b = zeros(2,2);
cps.cp_c = zeros(2,2);
cps.cp_d = zeros(2,2);
noElems = 16;
if nargin < 4
    elementList = 1:noElems;
end
index = [copyVector((1:4).',4,1), copyVector((1:4).',4,2)];
arr = pi*[0, 0.25, 0.5, 0.75, 1].'-pi/2;
arr(2) = etaShift-pi/2;
arr(end-1) = pi/2-etaShift;
elRangeEta = [arr(1:end-1), arr(2:end)];
arr = 2*pi*[0, 0.25, 0.5, 0.75, 1].'-pi;
elRangeXi = [arr(1:end-1), arr(2:end)];
for e = elementList
    if nargin < 5
        idXi = index(e,1);
        idEta = index(e,2);
        xiE = elRangeXi(idXi,:);
        etaE = elRangeEta(idEta,:);
    end
     
    maxItrs = 100;
    Eps = 1e4*eps;
    for i = 1:2
        x = xiE(i);
        dgdy = @(y) sphereMappingFarField(x, y, R, dmx, 'dgdy');
        dgdy2 = @(y) sphereMappingFarField(x, y, R, dmx, 'dgdy2');
        r = 0;
        if max(abs(dgdy(linspace(etaE(1),etaE(2),100)))) < Eps
            y_res = Inf;
            r = Inf;
        else
            y_res = newtonsMethod(dgdy,dgdy2,sum(etaE)/2,maxItrs,Eps);
        end
        if dgdy(y_res) > Eps
            y_res = Inf;
        end
         
        if ~isinf(y_res)
            if abs(dgdy2(y_res)) < eps
                r = 2;
            else
                r = 1;
            end
        end
        if i == 1
            cps(e).cp_a = [etaE(1), etaE(2);
                            0, 0];
            if etaE(1) < y_res && y_res < etaE(2)
                cps(e).cp_a = [etaE(1), y_res, etaE(2);
                               0, r, 0];
            elseif y_res == etaE(2)
                cps(e).cp_a(:,2) = [etaE(2); r];
            elseif y_res == etaE(1)
                cps(e).cp_a(:,1) = [etaE(1); r];
            elseif isinf(r)
                cps(e).cp_a = [etaE(1), etaE(2);
                                r, r];
            end
        else
            cps(e).cp_b = [etaE(1), etaE(2);
                            0, 0];
            if etaE(1) < y_res && y_res < etaE(2)
                cps(e).cp_b = [etaE(1), y_res, etaE(2);
                               0, r, 0];
            elseif y_res == etaE(1)
                cps(e).cp_b(:,1) = [etaE(1); r];
            elseif y_res == etaE(2)
                cps(e).cp_b(:,2) = [etaE(2); r];
            elseif isinf(r)
                cps(e).cp_b = [etaE(1), etaE(2);
                                r, r];
            end
        end
    end
    for i = 1:2
        y = etaE(i);
        dgdx = @(x) sphereMappingFarField(x, y, R, dmx, 'dgdx');
        dgdx2 = @(x) sphereMappingFarField(x, y, R, dmx, 'dgdx2');
        r = 0;
        if max(abs(dgdx(linspace(xiE(1),xiE(2),100)))) < Eps
            x_res = Inf;
            r = Inf;
        else
            x_res = newtonsMethod(dgdx,dgdx2,sum(xiE)/2,maxItrs,Eps);
        end
        if dgdx(x_res) > Eps
            x_res = Inf;
        end
         
        if ~isinf(x_res)
            if abs(dgdx2(x_res)) < eps
                r = 2;
            else
                r = 1;
            end
        end
        if i == 1
            cps(e).cp_c = [xiE(1), xiE(2);
                            0, 0];
            if xiE(1) < x_res && x_res < xiE(2)
                cps(e).cp_c = [xiE(1), x_res, xiE(2);
                               0, r, 0];
            elseif x_res == xiE(1)
                cps(e).cp_c(:,1) = [xiE(1); r];
            elseif x_res == xiE(2)
                cps(e).cp_c(:,2) = [xiE(2); r];
            elseif isinf(r)
                cps(e).cp_c = [xiE(1), xiE(2);
                                r, r];
            end
        else
            cps(e).cp_d = [xiE(1), xiE(2);
                            0, 0];
            if xiE(1) < x_res && x_res < xiE(2)
                cps(e).cp_d = [xiE(1), x_res, xiE(2);
                               0, r, 0];
            elseif x_res == xiE(1)
                cps(e).cp_d(:,1) = [xiE(1); r];
            elseif x_res == xiE(2)
                cps(e).cp_d(:,2) = [xiE(2); r];
            elseif isinf(r)
                cps(e).cp_d = [xiE(1), xiE(2);
                                r, r];
            end
        end
    end
end
if numel(elementList) == 1
    cps = cps(e);
end
cps.yx = true;
cps.xy = true;
if any(isinf(cps.cp_a(2,:))) || any(isinf(cps.cp_b(2,:)))
    cps.yx = false;
end
if any(isinf(cps.cp_c(2,:))) || any(isinf(cps.cp_d(2,:)))
    cps.xy = false;
end