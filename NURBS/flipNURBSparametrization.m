function nurbs = flipNURBSparametrization(nurbs,dir)
if iscell(nurbs)
    for patch = 1:numel(nurbs)
        coeffs = nurbs{patch}.coeffs;
        knots = nurbs{patch}.knots;
        switch nurbs{patch}.type
            case '3Dsurface'
                switch dir
                    case 'xi'
                        coeffs = coeffs(:,end:-1:1,:);
                        knots{1} = 1-knots{1}(end:-1:1);
                    case 'eta'
                        coeffs = coeffs(:,:,end:-1:1);
                        knots{2} = 1-knots{2}(end:-1:1);
                end
        end
        nurbs{patch}.coeffs = coeffs;
        nurbs{patch}.knots = knots;
    end
else
    switch dir
        case 'xi'
            nurbs.coeffs = nurbs.coeffs(:,end:-1:1,:);
            nurbs.knots{1} = 1-nurbs.knots{1}(end:-1:1);
        case 'eta'
            nurbs.coeffs = nurbs.coeffs(:,:,end:-1:1);
            nurbs.knots{2} = 1-nurbs.knots{2}(end:-1:1);
    end
end