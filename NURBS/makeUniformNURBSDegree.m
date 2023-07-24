function nurbs = makeUniformNURBSDegree(nurbs, p)

p_max = findMaxDegree(nurbs);
if nargin > 1
    p_max = max(p_max,p);
end
for i = 1:numel(nurbs)
    nurbs(i) = elevateNURBSdegree(nurbs(i),max(p_max-nurbs{i}.degree,0));
end