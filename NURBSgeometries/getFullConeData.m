function nurbs = getFullConeData(R_1, R_2, L, translateX)

nurbs1  = rotateNURBS(getConeData(R_1, R_2, L, 0, 2*pi/3, translateX),-pi/2,'Xaxis');
nurbs2 = rotateNURBS(nurbs1,2*pi/3,'Xaxis');
nurbs3 = rotateNURBS(nurbs2,2*pi/3,'Xaxis');
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
Eta = [0,0,1,1];
controlPts = zeros(4,7,2);
controlPts(:,1:3,:) = nurbs1{1}.coeffs;
controlPts(:,4:5,:) = nurbs2{1}.coeffs(:,2:3,:);
controlPts(:,6:7,:) = nurbs3{1}.coeffs(:,2:3,:);
nurbs = createNURBSobject(controlPts,{Xi, Eta});