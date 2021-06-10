function nurbs = explodeNURBS(nurbs,dir,reverseDir)

d_p = nurbs{1}.d_p;
if nargin < 2
    dir = 1:d_p;
end
if nargin < 3
    reverseDir = false;
end
if reverseDir
    j_dirs = numel(dir):-1:1;
else
    j_dirs = 1:numel(dir);
end
for j = j_dirs
    nurbs = explodeNURBSdir(nurbs,dir(j),d_p);
end

function nurbs2 = explodeNURBSdir(nurbs,dir,d_p)
nurbs2 = {};
for i = 1:numel(nurbs)
    nurbs_i = nurbs{i};
    knots = nurbs_i.knots;
    Xi = nurbs_i.knots{dir};
    p = nurbs_i.degree(dir);
    idx = 2;
    idxPrev = 2;
    while idx < numel(Xi)-p
        if Xi(idx+p) == Xi(idx+2*p-1)
            Xi2 = Xi(idxPrev:idx+2*p-1);
            xi_s = Xi2(1);
            xi_e = Xi2(end);
            Xi2 = [0, (Xi2-xi_s)/(xi_e-xi_s), 1];
            knots{dir} = Xi2;
            if d_p == 2
                switch dir
                    case 1
                        controlPts = nurbs_i.coeffs(:,idxPrev-1:idx+p-1,:);
                    case 2
                        controlPts = nurbs_i.coeffs(:,:,idxPrev-1:idx+p-1);
                end
            else
                switch dir
                    case 1
                        controlPts = nurbs_i.coeffs(:,idxPrev-1:idx+p-1,:,:);
                    case 2
                        controlPts = nurbs_i.coeffs(:,:,idxPrev-1:idx+p-1,:);
                    case 3
                        controlPts = nurbs_i.coeffs(:,:,:,idxPrev-1:idx+p-1);
                end
            end
            nurbs2{end+1} = nurbs_i;
            
            nurbs2{end}.coeffs = controlPts;
            nurbs2{end}.knots = knots;
            np = size(controlPts);
            nurbs2{end}.number = np(2:end);
    
            idx = idx+p-1;
            idxPrev = idx+1;
        end
        idx = idx + 1;
    end
end