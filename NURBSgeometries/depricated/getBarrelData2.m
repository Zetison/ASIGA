function nurbs = getBarrelData2(R_o, R_i, L)
error('Depricated use getBarrelData')

nurbs(1:3) = flipNURBSparametrization(getQuarterDiskData(R_o,'Xaxis',-[1,0,0]*L/2),'eta');
nurbs(4:6) = rotateNURBS(nurbs(1:3),pi/2,'Xaxis');
nurbs(7:9) = rotateNURBS(nurbs(1:3),2*pi/2,'Xaxis');
nurbs(10:12) = rotateNURBS(nurbs(1:3),3*pi/2,'Xaxis');

nurbs(21:23) = getQuarterDiskData(R_o,'Xaxis',[1,0,0]*L/2);
nurbs(24:26) = rotateNURBS(nurbs(21:23),pi/2,'Xaxis');
nurbs(27:29) = rotateNURBS(nurbs(21:23),2*pi/2,'Xaxis');
nurbs(30:32) = rotateNURBS(nurbs(21:23),3*pi/2,'Xaxis');

Xi = [0,0,0,1,1,1];
Eta = [0,0,1,1];
map1 = [2,3,5,6,8,9,11,12];
map2 = map1+20;
for i = 1:8
    controlPts = zeros(4,3,2);
    controlPts(:,:,1) = nurbs{map1(i)}.coeffs(:,:,end);
    controlPts(:,:,2) = nurbs{map2(i)}.coeffs(:,:,1);
    nurbs{i+12} = createNURBSobject(controlPts,{Xi, Eta});
end

nurbs(13:20) = elevateDegreeInPatches(nurbs(13:20),[0 1]);
nurbs(13:20) = insertKnotsInPatches(nurbs(13:20),0,3);


