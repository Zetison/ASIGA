function waterMesh = embedSolidInProlateSpheroid3(nurbs,nurbs2,c_z,R,alignWithAxis,x_0,alpha,take_nurbs2at_zeta1, ellipse3)

if nargin < 7
    alpha = 0;
end
if nargin < 8
    take_nurbs2at_zeta1 = true;
end
if nargin < 9
    ellipse3 = false;
end

Xi_solid = nurbs.knots{1};
Eta_solid = nurbs.knots{2};
Xi_solid2 = nurbs2.knots{1};
Eta_solid2 = nurbs2.knots{2};

if ellipse3
    ellipsoid = getEllipsoidalData3(R,c_z,alignWithAxis,x_0);
else
    ellipsoid = getEllipsoidalData(R,R,c_z,alignWithAxis,x_0);
end

Xi_ellipsoid = ellipsoid.knots{1};
Eta_ellipsoid = ellipsoid.knots{2};

Xi = mergeKnotVectors(Xi_solid, Xi_ellipsoid);
Eta = mergeKnotVectors(Eta_solid, Eta_ellipsoid);

newXiKnotsInEllipsoid = knotdiff(Xi,Xi_ellipsoid);
newEtaKnotsInEllipsoid = knotdiff(Eta,Eta_ellipsoid);

newXiKnotsInEllipsoid2 = knotdiff(Xi,Xi_solid2);
newEtaKnotsInEllipsoid2 = knotdiff(Eta,Eta_solid2);

newXiKnotsInSolid = knotdiff(Xi,Xi_solid);
newEtaKnotsInSolid = knotdiff(Eta,Eta_solid);


ellipsoid = insertKnotsInNURBS(ellipsoid,{newXiKnotsInEllipsoid newEtaKnotsInEllipsoid []});


R_x = rotationXaxis(alpha);

for i = 1:size(ellipsoid.coeffs,2)
    for j = 1:size(ellipsoid.coeffs,3)
        ellipsoid.coeffs(1:3,i,j) = R_x*ellipsoid.coeffs(1:3,i,j);
    end
end
         

nurbs = insertKnotsInNURBS(nurbs,{newXiKnotsInSolid newEtaKnotsInSolid []});
nurbs2 = insertKnotsInNURBS(nurbs2,{newXiKnotsInEllipsoid2 newEtaKnotsInEllipsoid2 []});

Zeta = [0 0 0 1 1 1];

controlPts = nurbs.coeffs(:,:,:,end);
if take_nurbs2at_zeta1
    controlPts(:,:,:,end+1) = nurbs2.coeffs(:,:,:,end);
else
    controlPts(:,:,:,end+1) = nurbs2.coeffs(:,:,:,1);
end
controlPts(:,:,:,end+1) = ellipsoid.coeffs;
          



waterMesh = createNURBSobject(controlPts,{Xi, Eta, Zeta});



function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];
