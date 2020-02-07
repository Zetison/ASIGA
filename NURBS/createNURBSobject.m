function nurbs = createNURBSobject(coeffs,knots)

np = size(coeffs);
switch length(knots)
    case 1
        if size(coeffs,1) == 2
            % constructing NURBS in 1D
            nurbs.type = '1Dnurbs';
            nurbs.coeffs = coeffs;
            n = np(2);
            nurbs.number = n;
            nurbs.degree = size(knots{1},2)-n-1;
            nurbs.knots = {knots{1}/knots{1}(end)};
        elseif size(coeffs,1) == 3
            % constructing a NURBS curve in 2D
            nurbs.type = '2Dcurve';
            nurbs.coeffs = coeffs;
            n = np(2);
            nurbs.number = n;
            nurbs.degree = size(knots{1},2)-n-1;
            nurbs.knots = {knots{1}/knots{1}(end)};
        elseif size(coeffs,1) == 4
            % constructing a NURBS curve in 3D
            nurbs.type = '3Dcurve';
            nurbs.coeffs = coeffs;
            n = np(2);
            nurbs.number = n;
            nurbs.degree = size(knots{1},2)-n-1;
            nurbs.knots = {knots{1}/knots{1}(end)};
        end
    case 2
        if size(coeffs,1) == 4
            % Construct NURBS surface in 3D
            nurbs.type = '3Dsurface';
            n = np(2);
            m = np(3);
            
            nurbs.coeffs = coeffs;
            
            Xi = knots{1};
            Eta = knots{2}; 
            nurbs.knots = {Xi/Xi(end), Eta/Eta(end)};
            
            p = size(knots{1},2)-n-1;
            q = size(knots{2},2)-m-1;
            nurbs.degree = [p q];
            nurbs.number = [n m];
        elseif size(coeffs,1) == 3
            % Construct NURBS surface in 2D
            nurbs.type = '2Dsurface';
            n = np(2);
            m = np(3);
            
            nurbs.coeffs = coeffs;
            
            Xi = knots{1};
            Eta = knots{2}; 
            nurbs.knots = {Xi/Xi(end), Eta/Eta(end)};
            
            p = size(knots{1},2)-n-1;
            q = size(knots{2},2)-m-1;
            nurbs.degree = [p q];
            nurbs.number = [n m];
        end
    case 3
        if size(coeffs,1) == 4
            % Construct NURBS volume in 3D
            nurbs.type = '3Dvolume';
            n = np(2);
            m = np(3);
            l = np(4);
            
            nurbs.coeffs = coeffs;
            
            Xi = knots{1};
            Eta = knots{2};
            Zeta = knots{3};    
            nurbs.knots = {Xi/Xi(end), Eta/Eta(end), Zeta/Zeta(end)}; % normalize knots
            p = size(knots{1},2)-n-1;
            q = size(knots{2},2)-m-1;
            r = size(knots{3},2)-l-1;
            nurbs.degree = [p q r];
            nurbs.number = [n m l];
        end
end
