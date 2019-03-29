function nurbs = getRodData2(leftStartingPoint,L,n,p)

Xi = [zeros(1,p+1) linspace2(0,1,n-(p+1)) ones(1,p+1)];

controlPts = zeros(2,n);
controlPts(1,:) = linspace(leftStartingPoint,leftStartingPoint+L,n);
controlPts(2,:) = ones(1,n);


nurbs.type = '1Dnurbs';
nurbs.coeffs = controlPts;
nurbs.number = n;
nurbs.degree = p;
nurbs.knots = Xi/Xi(end);