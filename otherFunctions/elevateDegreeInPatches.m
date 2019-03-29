function nurbs = elevateDegreeInPatches(nurbs,degreeElevArr)

if ~iscell(nurbs)
    nurbs = {nurbs};
end

for i = 1:numel(nurbs)  
    nurbs{i} = elevateNURBSdegree(nurbs{i},degreeElevArr);
end