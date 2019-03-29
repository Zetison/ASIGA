function solidAndEllipsoid = embedSolidInEllipse(solid,c_x,R)

Xi_solid = solid.knots;

ellipse = getEllipseData(c_x,R);

Xi_ellipse = ellipse.knots;

Xi = mergeKnotVectors(Xi_solid, Xi_ellipse);

newXiKnotsInEllipsoid = knotdiff(Xi,Xi_ellipse);

newXiKnotsInSolid = knotdiff(Xi,Xi_solid);


ellipse = insertKnotsInNURBS(ellipse,newXiKnotsInEllipsoid);

solid = insertKnotsInNURBS(solid,newXiKnotsInSolid);

Eta = [0 0 1 1];

controlPts = zeros(3,size(solid.coeffs,2), 2);
controlPts(:,:,1) = solid.coeffs;
controlPts(:,:,2) = ellipse.coeffs;
          


solidAndEllipsoid = createNURBSobject(controlPts,{Xi, Eta});

