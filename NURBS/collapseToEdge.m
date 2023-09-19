function nurbs = collapseToEdge(nurbs,midx,midx_midx)
% Collapses face midx of nurbs onto edge the masteridx midx_midx of face midx

coeffs = nurbs.coeffs;
knots = nurbs.knots;
degree = nurbs.degree;
switch midx
    case {1,2} % xi == 0 || xi == 1
        g = aveknt(knots{1},degree(1)+1);
        if midx == 2
            g = 1-g; % flip
        end
        switch midx_midx
            case 1 % Collapse to eta == 0
                coeffs = coeffs(:,:,1,:) + g.*(coeffs - coeffs(:,:,1,:));
            case 2 % Collapse to eta == 1
                coeffs = coeffs(:,:,end,:) + g.*(coeffs - coeffs(:,:,end,:));
            case 3 % Collapse to zeta == 0
                coeffs = coeffs(:,:,:,1) + g.*(coeffs - coeffs(:,:,:,1));
            case 4 % Collapse to zeta == 1
                coeffs = coeffs(:,:,:,end) + g.*(coeffs - coeffs(:,:,:,end));
        end
    case {3,4} % eta == 0 || eta == 1
        g = reshape(aveknt(knots{2},degree(2)+1),1,1,[]);
        if midx == 4
            g = 1-g; % flip
        end
        switch midx_midx
            case 1 % Collapse to zeta == 0
                coeffs = coeffs(:,:,:,1) + g.*(coeffs - coeffs(:,:,:,1));
            case 2 % Collapse to zeta == 1
                coeffs = coeffs(:,:,:,end) + g.*(coeffs - coeffs(:,:,:,end));
            case 3 % Collapse to xi == 0
                coeffs = coeffs(:,1,:,:) + g.*(coeffs - coeffs(:,1,:,:));
            case 4 % Collapse to xi == 1
                coeffs = coeffs(:,end,:,:) + g.*(coeffs - coeffs(:,end,:,:));
        end
    case {5,6} % zeta == 0 || zeta == 1
        g = reshape(aveknt(knots{3},degree(3)+1),1,1,1,[]);
        if midx == 6
            g = 1-g; % flip
        end
        switch midx_midx
            case 1 % Collapse to xi == 0
                coeffs = coeffs(:,1,:,:) + g.*(coeffs - coeffs(:,1,:,:));
            case 2 % Collapse to xi == 1
                coeffs = coeffs(:,end,:,:) + g.*(coeffs - coeffs(:,end,:,:));
            case 3 % Collapse to eta == 0
                coeffs = coeffs(:,:,1,:) + g.*(coeffs - coeffs(:,:,1,:));
            case 4 % Collapse to eta == 1
                coeffs = coeffs(:,:,end,:) + g.*(coeffs - coeffs(:,:,end,:));
        end
end
nurbs.coeffs = coeffs;