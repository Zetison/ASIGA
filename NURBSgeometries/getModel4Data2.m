function nurbs = getModel4Data2(R, t)

nurbs(1:3) = flipNURBSparametrization(getQuarterDiskData(R+t/2,'Xaxis',-[1,1,1]*t/2),'eta');
nurbs(4:6) = flipNURBSparametrization(getQuarterDiskData(R+t/2,'Yaxis',-[1,1,1]*t/2),'eta');
nurbs(7:9) = flipNURBSparametrization(getQuarterDiskData(R+t/2,'Zaxis',-[1,1,1]*t/2),'eta');

nurbs(16:18) = getQuarterDiskData(R-t/2,'Xaxis',[1,1,1]*t/2);
nurbs(19:21) = getQuarterDiskData(R-t/2,'Yaxis',[1,1,1]*t/2);
nurbs(22:24) = getQuarterDiskData(R-t/2,'Zaxis',[1,1,1]*t/2);

Xi = [0,0,0,1,1,1];
Eta = [0,0,1,1];
map1 = [2,3,5,6,8,9];
map2 = map1+15;
for i = 1:6
    controlPts = zeros(4,3,2);
    controlPts(:,:,1) = nurbs{map1(i)}.coeffs(:,:,end);
    controlPts(:,:,2) = nurbs{map2(i)}.coeffs(:,:,1);
    nurbs{i+9} = createNURBSobject(controlPts,{Xi, Eta});
end

nurbs(10:15) = elevateDegreeInPatches(nurbs(10:15),[0 1]);


