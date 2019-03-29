function nurbs = flipNURBSparametrization(nurbs,dir)

switch dir
    case 'xi'
        nurbs.coeffs = nurbs.coeffs(:,end:-1:1,:);
        nurbs.knots{1} = 1-nurbs.knots{1}(end:-1:1);
    case 'eta'
        nurbs.coeffs = nurbs.coeffs(:,:,end:-1:1);
        nurbs.knots{2} = 1-nurbs.knots{2}(end:-1:1);
end