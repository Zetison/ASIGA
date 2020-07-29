function nurbs = permuteNURBS(nurbs,order)

% rearranges the dimensions of the parametric direction in nurbs so that
% they are in the order specified by the vector order


if ~iscell(nurbs)
    nurbs = {nurbs};
end

for patch = 1:numel(nurbs)
    nurbs{patch}.coeffs = permute(nurbs{patch}.coeffs,[1,order+1]);
    nurbs{patch}.knots = nurbs{patch}.knots(order);
    nurbs{patch}.number = nurbs{patch}.number(order);
    nurbs{patch}.degree = nurbs{patch}.degree(order);
end