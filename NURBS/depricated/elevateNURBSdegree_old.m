function newNurbs = elevateNURBSdegree(nurbs,degreeElevations)
error('Depricated')
if strcmp(nurbs.type, '3Dvolume')
    % NURBS represents a volume
    n = nurbs.number(1);
    m = nurbs.number(2);
    l = nurbs.number(3);

    coeffs = nurbs.coeffs;
    for i = 1:n
        for j = 1:m
            for k = 1:l
                coeffs(1:3,i,j,k) = coeffs(1:3,i,j,k)*coeffs(4,i,j,k);
            end
        end
    end
    % Elevate in degree in xi direction
    if degreeElevations(1) == 0;
        knots{1} = nurbs.knots{1};
    else   
        coeffs = permute(coeffs,[1 3 4 2]);
        coeffs = reshape(coeffs,4*m*l,n);
        [coeffs, knots{1}] = elevateBsplinesDegree(nurbs.number(1), nurbs.degree(1), nurbs.knots{1}, coeffs, degreeElevations(1));
        n = size(coeffs,2);
        coeffs = reshape(coeffs,[4 m l n]);
        coeffs = permute(coeffs,[1 4 2 3]);
    end
    % Elevate in degree in eta direction
    if degreeElevations(2) == 0;
        knots{2} = nurbs.knots{2};
    else   
        coeffs = permute(coeffs,[1 2 4 3]);
        coeffs = reshape(coeffs,4*n*l,m);
        [coeffs, knots{2}] = elevateBsplinesDegree(nurbs.number(2), nurbs.degree(2), nurbs.knots{2}, coeffs, degreeElevations(2));
        m = size(coeffs,2);
        coeffs = reshape(coeffs,[4 n l m]);
        coeffs = permute(coeffs,[1 2 4 3]);
    end
    % Elevate in degree in zeta direction
    if degreeElevations(3) == 0;
        knots{3} = nurbs.knots{3};
    else   
        coeffs = reshape(coeffs,4*n*m,l);
        [coeffs, knots{3}] = elevateBsplinesDegree(nurbs.number(3), nurbs.degree(3), nurbs.knots{3}, coeffs, degreeElevations(3));
        l = size(coeffs,2);
        coeffs = reshape(coeffs,[4 n m l]);
    end
    % Finale projection step
    for i = 1:n
        for j = 1:m
            for k = 1:l
                coeffs(1:3,i,j,k) = coeffs(1:3,i,j,k)/coeffs(4,i,j,k);
            end
        end
    end
elseif strcmp(nurbs.type, '2Dcurve')
    % NURBS represents a curve in 2D
    if degreeElevations == 0
        coeffs = nurbs.coefs;
        knots = nurbs.knots;
    else
        n = nurbs.number;
        coeffs = nurbs.coeffs';
        coeffs(:,1:2) = coeffs(:,1:2).*[coeffs(:,3) coeffs(:,3)];
        
        coeffs = reshape(coeffs,3,n);
        [coeffs,knots] = elevateBsplinesDegree(nurbs.number, nurbs.degree, nurbs.knots, coeffs, degreeElevations);
        n = size(coeffs,2);
        coeffs = reshape(coeffs,n,3);
        
        coeffs(:,1:2) = coeffs(:,1:2)./[coeffs(:,3) coeffs(:,3)];
        coeffs = coeffs';
    end
elseif strcmp(nurbs.type, '1Dnurbs') 
    % NURBS in 1D
    if degreeElevations == 0
        coeffs = nurbs.coefs;
        knots = nurbs.knots;
    else
        n = nurbs.number;
        coeffs = nurbs.coeffs';
        coeffs(:,1) = coeffs(:,1).*coeffs(:,2);
        
        coeffs = reshape(coeffs,2,n);
        [coeffs,knots] = elevateBsplinesDegree(nurbs.number, nurbs.degree, nurbs.knots, coeffs, degreeElevations);
        n = size(coeffs,2);
        coeffs = reshape(coeffs,n,2);
        
        coeffs(:,1) = coeffs(:,1)./coeffs(:,2);
        coeffs = coeffs';
    end
end

% construct new NURBS
newNurbs = createNURBSobject(coeffs,knots);










