function refinednurbs = refineNURBSevenly(nurbs,n)

refinednurbs = nurbs;
for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    for dir = 1:d_p
        uniqueXi = unique(nurbs{patch}.knots{dir});
        dir2 = setdiff(1:d_p,dir);
        uniqueEta = unique(nurbs{patch}.knots{dir2});
        parm_pts = NaN(numel(uniqueEta),d_p);
        parm_pts(:,dir2) = uniqueEta;
        for j = 1:numel(uniqueXi)-1
            I_max = -Inf;
            for i = 1:size(parm_pts,1)
                I = NURBSarcLength(nurbs{patch},uniqueXi(j),uniqueXi(j+1),parm_pts(i,:),dir);
                if I_max < I
                    I_max = I;
                end
            end
            newXiKnots = linspace2(uniqueXi(j),uniqueXi(j+1),round(I_max*n));
            newKnots = cell(1,d_p);
            newKnots{dir} = newXiKnots;
            refinednurbs(patch) = insertKnotsInNURBS(refinednurbs(patch),newKnots);
        end
    end
end