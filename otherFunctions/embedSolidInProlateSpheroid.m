function waterMesh = embedSolidInEllipsoid(nurbs,c_z,R,alignWithAxis,x_0, alpha)

if nargin < 6
    alpha = 0;
end

p_xi = nurbs.degree(1);
p_eta = nurbs.degree(2);
Xi_solid = nurbs.knots{1};
Eta_solid = nurbs.knots{2};

ellipsoid = getEllipsoidalData(a,b,c,alignWithAxis,x_0, alpha, Xi_solid(p_xi+2:end-p_xi-1), Eta_solid(p_eta+2:end-p_eta-1));

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

Zeta = [0 0 1 1];

R_x = rotationXaxis(alpha);

controlPts = nurbs.coeffs(:,:,:,end);

for i = 1:size(ellipsoid.coeffs,2)
    for j = 1:size(ellipsoid.coeffs,3)
        ellipsoid.coeffs(1:3,i,j) = R_x*ellipsoid.coeffs(1:3,i,j);
    end
end
controlPts(:,:,:,end+1) = ellipsoid.coeffs;
          



waterMesh = createNURBSobject(controlPts,{Xi, Eta, Zeta});


function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];
