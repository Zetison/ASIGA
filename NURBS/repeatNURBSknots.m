function nurbs = repeatNURBSknots(nurbs)
if ~iscell(nurbs)
    nurbs = {nurbs};
end

for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    knots = nurbs{patch}.knots;
    degree = nurbs{patch}.degree;
    newKnots = cell(1,d_p);
    Eps = 1e5*eps;
    for j = 1:d_p  
        uniqueKnots = unique(knots{j});
        for i = 1:length(uniqueKnots)
            mm = length(find(abs(knots{j} - uniqueKnots(i)) < Eps));  
            newKnots{j} = [newKnots{j}, uniqueKnots(i)*ones(1,degree(j)-mm)];
        end
    end
    nurbs(patch) = insertKnotsInNURBS(nurbs(patch),newKnots);
end