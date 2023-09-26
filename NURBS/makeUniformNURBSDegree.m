function [nurbs,degreeElevations] = makeUniformNURBSDegree(nurbs, p)

p_max = findMaxDegree(nurbs);
if nargin > 1
    p_max = max(p_max,p);
end
degreeElevations = zeros(numel(nurbs),numel(p_max));
for i = 1:numel(nurbs)
    degreeElevations(i,:) = max(p_max-nurbs{i}.degree,0);
    nurbs(i) = elevateNURBSdegree(nurbs(i),degreeElevations);
end