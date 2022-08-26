function [hmax, hmin, diagsMax] = findMaxElementDiameter(patches)

if ~iscell(patches)
    patches = {patches};
end
H_max = zeros(numel(patches),1);
H_min = zeros(numel(patches),1);
diagsMax = cell(numel(patches),1);
for patch = 1:numel(patches)
    nurbs = patches{patch}.nurbs;
    d_p = nurbs.d_p;
    d = nurbs.d;
    uniqueKnots = cell(1,d_p);
    noKnots = zeros(1,d_p);
    for i = 1:d_p
        uniqueKnots{i} = unique(nurbs.knots{i});
        noKnots(i) = numel(uniqueKnots{i});
    end

    switch d_p
        case 1
            v = evaluateNURBSvec(nurbs,uniqueKnots{1}.');
            d_1 = vecnorm(v(2:end,:) - v(1:end-1,:), 2, d_p+1);
            diaMax = d_1;
            diaMin = d_1;
        case 2
            [XI,ETA] = ndgrid(uniqueKnots{1},uniqueKnots{2});

            v = evaluateNURBSvec(nurbs,[XI(:),ETA(:)]);
            v = reshape(v,[noKnots,d]);
            d_1 = vecnorm(v(2:end,2:end,:) - v(1:end-1,1:end-1,:), 2, d_p+1);
            d_2 = vecnorm(v(1:end-1,2:end,:) - v(2:end,1:end-1,:), 2, d_p+1);
            diaMax = max([d_1(:), d_2(:)],[],2);
            diaMin = min([d_1(:), d_2(:)],[],2);
        case 3
            [XI,ETA,ZETA] = ndgrid(uniqueKnots{1},uniqueKnots{2},uniqueKnots{3});

            v = evaluateNURBSvec(nurbs,[XI(:),ETA(:),ZETA(:)]);
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