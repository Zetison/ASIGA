function waterMesh = embedSolidInProlateSpheroid2(nurbs,nurbs2,nurbs3,nurbs4,c_z,R,alignWithAxis)

Xi_solid = nurbs.knots{1};
Eta_solid = nurbs.knots{2};

ellipsoid = getEllipsoidalData(R,R,c_z,alignWithAxis);


Xi_ellipsoid = ellipsoid.knots{1};
Eta_ellipsoid = ellipsoid.knots{2};

Xi = mergeKnotVectors(Xi_solid, Xi_ellipsoid);
Eta = mergeKnotVectors(Eta_solid, Eta_ellipsoid);

newXiKnotsInEllipsoid = knotdiff(Xi,Xi_ellipsoid);
newEtaKnotsInEllipsoid = knotdiff(Eta,Eta_ellipsoid);

newXiKnotsInSolid = knotdiff(Xi,Xi_solid);
newEtaKnotsInSolid = knotdiff(Eta,Eta_solid);


ellipsoid = insertKnotsInNURBS(ellipsoid,{newXiKnotsInEllipsoid newEtaKnotsInEllipsoid []});

nurbs = insertKnotsInNURBS(nurbs,{newXiKnotsInSolid newEtaKnotsInSolid []});
nurbs2 = insertKnotsInNURBS(nurbs2,{newXiKnotsInSolid newEtaKnotsInSolid []});
nurbs3 = insertKnotsInNURBS(nurbs3,{newXiKnotsInSolid newEtaKnotsInSolid []});
nurbs4 = insertKnotsInNURBS(nurbs4,{newXiKnotsInSolid newEtaKnotsInSolid []});

Zeta = [0 0 0 0.5 0.5 1 1 1];


controlPts = nurbs.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = nurbs2.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = nurbs3.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = nurbs4.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = ellipsoid.coeffs;
          



waterMesh = createNURBSobject(controlPts,{Xi, Eta, Zeta});

