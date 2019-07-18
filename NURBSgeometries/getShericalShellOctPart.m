function nurbs = getShericalShellOctPart(R)

rotAxis = [1,1,1]/sqrt(3);
nurbs{1} = getSphericalShellDataPatchedQuarter(R);
nurbs(2) = rotateNURBS(nurbs(1),2*pi/3, rotAxis);
nurbs(3) = rotateNURBS(nurbs(1),4*pi/3, rotAxis);