function [hmax, hmin, diagsMax] = findMaxElementDiameter(nurbs,useGrevPts)

if ~iscell(nurbs)
    nurbs = {nurbs};
end
if nargin < 2
    useGrevPts = false;
end
H_max = zeros(numel(nurbs),1);
H_min = zeros(numel(nurbs),1);
diagsMax = cell(numel(nurbs),1);
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    d = nurbs{patch}.d;
    uniqueKnots = cell(1,d_p);
    noKnots = zeros(1,d_p);
    for i = 1:d_p
        if useGrevPts
            uniqueKnots{i} = aveknt(nurbs{patch}.knots{i},nurbs{patch}.degree(i)+1);
        else
            uniqueKnots{i} = unique(nurbs{patch}.knots{i});
        end
        noKnots(i) = numel(uniqueKnots{i});
    end

    switch d_p
        case 1
            v = evaluateNURBSvec(nurbs{patch},uniqueKnots{1}.');
            d_1 = vecnorm(v(2:end,:) - v(1:end-1,:), 2, d_p+1);
            diaMax = d_1;
            diaMin = d_1;
        case 2
            [XI,ETA] = ndgrid(uniqueKnots{1},uniqueKnots{2});

            v = evaluateNURBSvec(nurbs{patch},[XI(:),ETA(:)]);
            v = reshape(v,[noKnots,d]);
            d_1 = vecnorm(v(2:end,2:end,:) - v(1:end-1,1:end-1,:), 2, d_p+1);
            d_2 = vecnorm(v(1:end-1,2:end,:) - v(2:end,1:end-1,:), 2, d_p+1);
            diaMax = max([d_1(:), d_2(:)],[],2);
            diaMin = min([d_1(:), d_2(:)],[],2);
        case 3
            [XI,ETA,ZETA] = ndgrid(uniqueKnots{1},uniqueKnots{2},uniqueKnots{3});

            v = evaluateNURBSvec(nurbs{patch},[XI(:),ETA(:),ZETA(:)]);
            v = reshape(v,[noKnots,d]);
            d_1 = vecnorm(v(2:end,2:end,2:end,:) - v(1:end-1,1:end-1,1:end-1,:), 2, d_p+1);
            d_2 = vecnorm(v(1:end-1,2:end,2:end,:) - v(2:end,1:end-1,1:end-1,:), 2, d_p+1);
            d_3 = vecnorm(v(2:end,1:end-1,2:end,:) - v(1:end-1,2:end,1:end-1,:), 2, d_p+1);
            d_4 = vecnorm(v(2:end,2:end,1:end-1,:) - v(1:end-1,1:end-1,2:end,:), 2, d_p+1);
            diaMax = max([d_1(:), d_2(:), d_3(:), d_4(:)],[],2);
            diaMin = min([d_1(:), d_2(:), d_3(:), d_4(:)],[],2);
    end
    H_max(patch) = max(diaMax);
    H_min(patch) = min(diaMin);
    diagsMax{patch} = diaMax;
end
diagsMax = cell2mat(diagsMax);
hmin = min(H_min);
hmax = max(H_max);