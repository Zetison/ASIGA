function nurbs = getQuarterAnnulus(R_i, R_o)
error('Depricated: use getDiskData')

Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
controlPts = zeros(3,3,2);

controlPts(:,1,1) = [ R_i  0       1           ];
controlPts(:,2,1) = [ R_i  R_i     1/sqrt(2)   ];
controlPts(:,3,1) = [ 0    R_i     1           ];

controlPts(:,1,2) = [ R_o  0       1           ];
controlPts(:,2,2) = [ R_o  R_o     1/sqrt(2)   ];
controlPts(:,3,2) = [ 0    R_o     1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});


