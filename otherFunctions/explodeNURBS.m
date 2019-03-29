function nurbsCol2 = explodeNURBS(nurbsCol,dir)
if ~iscell(nurbsCol)
    nurbsCol = {nurbsCol};
end
nurbsCol2 = {};
switch dir
    case 'xi'
        parm = 1;
    case 'eta'
        parm = 2;
    case 'zeta'
        parm = 3;
end
for i = 1:numel(nurbsCol)
    nurbs = nurbsCol{i};
    noParams = numel(nurbs.degree);
    knots = nurbs.knots;
    Xi = nurbs.knots{parm};
    p = nurbs.degree(parm);
    idx = 2;
    idxPrev = 2;
    while idx < numel(Xi)-p
        if Xi(idx+p) == Xi(idx+2*p-1)
            Xi2 = Xi(idxPrev:idx+2*p-1);
            xi_s = Xi2(1);
            xi_e = Xi2(end);
            Xi2 = [0, (Xi2-xi_s)/(xi_e-xi_s), 1];
            knots{parm} = Xi2;
            if noParams == 2
                if parm == 1
                    controlPts = nurbs.coeffs(:,idxPrev-1:idx+p-1,:);
                elseif parm == 2
                    controlPts = nurbs.coeffs(:,:,idxPrev-1:idx+p-1);
                end
            else
                if parm == 1
                    controlPts = nurbs.coeffs(:,idxPrev-1:idx+p-1,:,:);
                elseif parm == 2
                    controlPts = nurbs.coeffs(:,:,idxPrev-1:idx+p-1,:);
                elseif parm == 3
                    controlPts = nurbs.coeffs(:,:,:,idxPrev-1:idx+p-1);
                end
            end
            nurbsCol2{end+1} = createNURBSobject(controlPts,knots);
            idx = idx+p-1;
            idxPrev = idx+1;
        end
        idx = idx + 1;
    end
end