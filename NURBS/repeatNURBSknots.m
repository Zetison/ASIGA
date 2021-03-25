function nurbs = repeatNURBSknots(nurbs)
if ~iscell(nurbs)
    nurbs = {nurbs};
end

for patch = 1:numel(nurbs)
    d_p = nurbs{patch}.d_p;
    knots = nurbs{patch}.knots;
    degree = nurbs{patch}.degree;
    newKnots = cell(1,d_p);
    for j = 1:d_p  
        for i = 1:length(knots{j})
            mm = length(find(knots{j} == knots{j}(i)));  
            newKnots{j} = [newKnots{j}, knots{j}(i)*ones(degree(j)-mm,1)];
        end
    end
    nurbs(patch) = insertKnotsInNURBS(nurbs(patch),newKnots);
end