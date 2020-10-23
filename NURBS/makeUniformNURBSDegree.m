function nurbs = makeUniformNURBSDegree(nurbs, p)

d_p = numel(nurbs{1}.knots);
p_max = -Inf*ones(1,d_p);
for i = 1:numel(nurbs)
    for j = 1:d_p
        if nurbs{i}.degree(j) > p_max(j)
            p_max(j) = nurbs{i}.degree(j);
        end
    end
end
if nargin > 1
    p_max = max(p_max,p);
end
for i = 1:numel(nurbs)
    nurbs(i) = elevateNURBSdegree(nurbs(i),max(p_max-nurbs{i}.degree,0));
end