function nurbs = homogenizeNURBSparametrization(nurbs)


% Find maximal polynomial order and the union of all knot vectors
noPatches = numel(nurbs);
d_p = nurbs{1}.d_p; % Assume all patches to have the same parametric dimension
p_max = zeros(1,d_p);
unionKnots = cell(1,d_p);
for patch = 1:noPatches
    p_max = max(nurbs{patch}.degree,p_max);
    % Find the union of all knots
    for i = 1:d_p
        unionKnots{i} = [unionKnots{i}, setdiffUnique(unionKnots{i},nurbs{patch}.knots{i})];
    end
end
% Homogenize all nurbs
for patch = 1:noPatches
    newKnots = cell(1,d_p);
    for i = 1:d_p
        newKnots{i} = setdiffUnique(unionKnots{i},nurbs{patch}.knots{i});
    end
    nurbs(patch) = elevateNURBSdegree(nurbs(patch),p_max-nurbs{patch}.degree);
    nurbs(patch) = insertKnotsInNURBS(nurbs(patch),newKnots);
end

