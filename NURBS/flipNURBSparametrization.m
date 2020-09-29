function nurbs = flipNURBSparametrization(nurbs,dir)
if nargin < 2
    dir = 1;
end
for patch = 1:numel(nurbs)
    for i = 1:numel(dir)
        nurbs{patch}.coeffs = flip(nurbs{patch}.coeffs,dir(i)+1);
        nurbs{patch}.knots{dir(i)} = 1-nurbs{patch}.knots{dir(i)}(end:-1:1);
    end
end