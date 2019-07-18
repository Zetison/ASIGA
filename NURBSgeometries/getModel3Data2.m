function nurbs = getModel3Data2(R_1, R_2, L)

rotAxis = 'Xaxis';
nurbs = getShericalShellOctPart(R_1);
nurbs(4:6) = rotateNURBS(nurbs(1:3),pi/2, rotAxis);
nurbs(7:9) = rotateNURBS(nurbs(1:3),2*pi/2, rotAxis);
nurbs(10:12) = rotateNURBS(nurbs(1:3),3*pi/2, rotAxis);
nurbs(1:12) = translateNURBS(mirrorNURBS(nurbs,'x'),[-L,0,0]);

nurbs(21:23) = getShericalShellOctPart(R_2);
nurbs(24:26) = rotateNURBS(nurbs(21:23),pi/2, rotAxis);
nurbs(27:29) = rotateNURBS(nurbs(21:23),2*pi/2, rotAxis);
nurbs(30:32) = rotateNURBS(nurbs(21:23),3*pi/2, rotAxis);


Xi = [0,0,0,0,0,1,1,1,1,1];
Eta = [0,0,1,1];
map1 = [1,3,4,6,7,9,10,12];
map2 = map1+20;
for i = 1:2:8
    controlPts = zeros(4,5,2);
    controlPts(:,:,1) = nurbs{map1(i)}.coeffs(:,1,:);
    controlPts(:,:,2) = nurbs{map2(i)}.coeffs(:,1,:);
    nurbs{i+12} = createNURBSobject(controlPts,{Xi, Eta});
end
for i = 2:2:8
    controlPts = zeros(4,5,2);
    controlPts(:,:,1) = nurbs{map1(i)}.coeffs(:,:,1);
    controlPts(:,:,2) = nurbs{map2(i)}.coeffs(:,:,1);
    nurbs{i+12} = createNURBSobject(controlPts,{Xi, Eta});
end

nurbs(13:20) = elevateDegreeInPatches(nurbs(13:20),[0 3 0]);

