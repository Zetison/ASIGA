function nurbs = getCylinderData(R, L)

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 1 1];

controlPts = zeros(4,9,2);

controlPts(:,1,1,1) = [  R	0       0	 1           ];
controlPts(:,2,1,1) = [  R	R     0	 1/sqrt(2)   ];
controlPts(:,3,1,1) = [  0      R     0    1           ];
controlPts(:,4,1,1) = [ -R	R     0    1/sqrt(2)   ];
controlPts(:,5,1,1) = [ -R    0       0    1           ];
controlPts(:,6,1,1) = [ -R   -R     0    1/sqrt(2)   ];
controlPts(:,7,1,1) = [  0     -R     0    1           ];
controlPts(:,8,1,1) = [  R   -R     0    1/sqrt(2)   ];
controlPts(:,9,1,1) = [  R	0       0    1           ];

controlPts(:,1,2,1) = [  R	0       L    1           ];
controlPts(:,2,2,1) = [  R	R     L    1/sqrt(2)   ];
controlPts(:,3,2,1) = [  0  	R     L    1           ];
controlPts(:,4,2,1) = [ -R    R     L    1/sqrt(2)   ];
controlPts(:,5,2,1) = [ -R    0       L    1           ];
controlPts(:,6,2,1) = [ -R   -R     L    1/sqrt(2)   ];
controlPts(:,7,2,1) = [  0     -R     L    1           ];
controlPts(:,8,2,1) = [  R   -R     L    1/sqrt(2)   ];
controlPts(:,9,2,1) = [  R  	0       L    1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});