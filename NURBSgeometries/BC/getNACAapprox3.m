function cntrlPts = getNACAapprox3(t,p,Xi,f,weigh_ss,idx,xValues,xIndices, axisFlipIdx, axisFlipIdx2)

% f = @(xi) [xi.^2, getNACA(xi,t)];
enforceFrontSmoothness = true; % This was set to false in the original BCA approximations

if nargin < 6 || isempty(xValues)
    [xi_rep, I] = getRepeatedKnots(Xi,p);
    if enforceFrontSmoothness
        dofsToRemove{1} = sort(unique([2,I]));
    else
        dofsToRemove{1} = I;
    end
    dofsToRemove{2} = I;
    dofsToRemove{3} = I;
    temp = f(xi_rep.');
    if enforceFrontSmoothness
        values{1} = [temp(1,1); temp(:,1)];
    else
        values{1} = temp(:,1);
    end
    values{2} = temp(:,2);
    values{3} = temp(:,3);
    nurbs = leastSquares1D(Xi,p,f,3,weigh_ss,dofsToRemove,values);
else
    [xi_rep, I] = getRepeatedKnots(Xi,p);
    if enforceFrontSmoothness
        dofsToRemove{1} = sort(unique([2,I,xIndices]));
    else
        dofsToRemove{1} = sort(unique([I,xIndices]));
    end
    dofsToRemove{axisFlipIdx} = sort([I,I(idx)-1,I(idx)+1]);
    dofsToRemove{axisFlipIdx2} = I;
    temp = f(xi_rep.');
    values{1} = temp(:,1);
    values{2} = temp(:,2);
    values{3} = temp(:,3);
    values{axisFlipIdx} = [temp(1:idx-1,axisFlipIdx); temp(idx,axisFlipIdx)*ones(3,1); temp(idx+1:end,axisFlipIdx)];
    if enforceFrontSmoothness
        values{1} = [temp(1,1); values{1}(1:idx-1); xValues; values{1}(idx+1:end)];
    else
        values{1} = [values{1}(1:idx-1); xValues; values{1}(idx+1:end)];
    end
    nurbs = leastSquares1D(Xi,p,f,3,weigh_ss,dofsToRemove,values);
end
cntrlPts = nurbs{1}.coeffs(1:3,:).';