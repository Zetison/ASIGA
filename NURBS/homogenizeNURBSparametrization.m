function [nurbs,newKnots,degreeElevations] = homogenizeNURBSparametrization(nurbs)
% Ensure that the basis is the same for all patches by order elevation and knot insertion

p_max = findMaxDegree(nurbs);
[nurbs,degreeElevations] = makeUniformNURBSDegree(nurbs,p_max);
noPatches = numel(nurbs);

% Find the union of all knots
knots = nurbs{1}.knots;
for i = 2:noPatches
    for j = 1:nurbs{i}.d_p
        knots{j} = sort([knots{j}, setdiffUnique(nurbs{i}.knots{j},knots{j})]);
    end
end
% Homogenize all nurbs
newKnots = cell(1,noPatches);
for i = 1:noPatches
    newKnots{i} = cell(1,nurbs{i}.d_p);
    for j = 1:nurbs{i}.d_p
        newKnots{i}{j} = setdiffUnique(knots{j},nurbs{i}.knots{j});
    end
    nurbs(i) = insertKnotsInNURBS(nurbs(i),newKnots{i});
end