function nurbs = makeUniformNURBSDegree(nurbs, p)

p_max = -Inf;
for i = 1:numel(nurbs)
    if max(nurbs{i}.degree) > p_max
        p_max = max(nurbs{i}.degree);
    end
end
if nargin > 1
    p_max = max(p_max,p);
end
for i = 1:numel(nurbs)
    d_p = numel(nurbs{i}.knots);
    nurbs(i) = elevateNURBSdegree(nurbs(i),max(p_max*ones(1,d_p)-nurbs{i}.degree,0));
end