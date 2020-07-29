function nurbs = flipNURBSparametrization(nurbs,dir)
if nargin < 2
    dir = 1;
end
for patch = 1:numel(nurbs)
    nurbs{patch}.coeffs = flip(nurbs{patch}.coeffs,dir+1);
    nurbs{patch}.knots{dir} = 1-nurbs{patch}.knots{dir}(end:-1:1);
end