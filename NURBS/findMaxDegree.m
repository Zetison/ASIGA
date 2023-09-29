function degree = findMaxDegree(nurbs)

d_p = numel(nurbs{1}.knots);
degree = -Inf*ones(1,d_p);
for i = 1:numel(nurbs)
    for j = 1:d_p
        if nurbs{i}.degree(j) > degree(j)
            degree(j) = nurbs{i}.degree(j);
        end
    end
end