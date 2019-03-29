function nurbs = getRectangleIn3DData(P_1,P_2,P_3,P_4)


Xi = [0 0 1 1];
Eta = [0 0 1 1];

controlPts = ones(4,2,2);


controlPts(1:3,1,1) = P_1;
controlPts(1:3,2,1) = P_2;

controlPts(1:3,1,2) = P_4;
controlPts(1:3,2,2) = P_3;

nurbs = createNURBSobject(controlPts,{Xi, Eta});
