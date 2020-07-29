function nurbsSurf = extractSurface(nurbs, dir, surface)
error('Depricated. Use subNURBS() instead')

if nargin < 3
    dir = 'zeta';
end
if ~iscell(nurbs)
    nurbs = {nurbs};
end
switch surface
    case 'inner'
        switch dir
            case 'zeta'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Xi = nurbs{i}.knots{1};
                    Eta = nurbs{i}.knots{2};

                    controlPts = nurbs{i}.coeffs(:,:,:,1);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Xi, Eta});
                end
            case 'eta'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Xi = nurbs{i}.knots{1};
                    Zeta = nurbs{i}.knots{3};
                    n_xi = nurbs{i}.number(1);
                    n_zeta = nurbs{i}.number(3);
                    controlPts = reshape(nurbs{i}.coeffs(:,:,1,:),4,n_xi,n_zeta);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Xi, Zeta});
                end
            case 'xi'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Eta = nurbs{i}.knots{2};
                    Zeta = nurbs{i}.knots{3};
                    n_eta = nurbs{i}.number(2);
                    n_zeta = nurbs{i}.number(3);
                    controlPts = reshape(nurbs{i}.coeffs(:,1,:,:),4,n_eta,n_zeta);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Eta, Zeta});
                end
        end
    case 'outer'
        switch dir
            case 'zeta'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Xi = nurbs{i}.knots{1};
                    Eta = nurbs{i}.knots{2};

                    controlPts = nurbs{i}.coeffs(:,:,:,end);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Xi, Eta});
                end
            case 'eta'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Xi = nurbs{i}.knots{1};
                    Zeta = nurbs{i}.knots{3};
                    n_xi = nurbs{i}.number(1);
                    n_zeta = nurbs{i}.number(3);
                    controlPts = reshape(nurbs{i}.coeffs(:,:,end,:),4,n_xi,n_zeta);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Xi, Zeta});
                end
            case 'xi'
                nurbsSurf = cell(1,numel(nurbs));
                for i = 1:numel(nurbs)
                    Eta = nurbs{i}.knots{2};
                    Zeta = nurbs{i}.knots{3};
                    n_eta = nurbs{i}.number(2);
                    n_zeta = nurbs{i}.number(3);
                    controlPts = reshape(nurbs{i}.coeffs(:,end,:,:),4,n_eta,n_zeta);

                    nurbsSurf{i} = createNURBSobject(controlPts,{Eta, Zeta});
                end
        end
end