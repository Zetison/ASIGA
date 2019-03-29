function solid = getCircularPlateData(R, t)

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,9,2,2);

controlPts(:,1,1,1) = [  0      0       0    1           ];
controlPts(:,2,1,1) = [  0      0       0    1/sqrt(2)   ];
controlPts(:,3,1,1) = [  0      0       0    1           ];
controlPts(:,4,1,1) = [  0      0       0    1/sqrt(2)   ];
controlPts(:,5,1,1) = [  0      0       0    1           ];
controlPts(:,6,1,1) = [  0      0       0    1/sqrt(2)   ];
controlPts(:,7,1,1) = [  0      0       0    1           ];
controlPts(:,8,1,1) = [  0      0       0    1/sqrt(2)   ];
controlPts(:,9,1,1) = [  0      0       0    1           ];

controlPts(:,1,2,1) = [  R      0       0    1           ];
controlPts(:,2,2,1) = [  R      R       0    1/sqrt(2)   ];
controlPts(:,3,2,1) = [  0   	R       0    1           ];
controlPts(:,4,2,1) = [ -R      R       0    1/sqrt(2)   ];
controlPts(:,5,2,1) = [ -R      0       0    1           ];
controlPts(:,6,2,1) = [ -R     -R       0    1/sqrt(2)   ];
controlPts(:,7,2,1) = [  0     -R       0    1           ];
controlPts(:,8,2,1) = [  R     -R       0    1/sqrt(2)   ];
controlPts(:,9,2,1) = [  R      0       0    1           ];

controlPts(:,1,1,2) = [  0      0       t    1           ];
controlPts(:,2,1,2) = [  0      0       t    1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0      0       t    1           ];
controlPts(:,4,1,2) = [  0      0       t    1/sqrt(2)   ];
controlPts(:,5,1,2) = [  0      0       t    1           ];
controlPts(:,6,1,2) = [  0      0       t    1/sqrt(2)   ];
controlPts(:,7,1,2) = [  0      0       t    1           ];
controlPts(:,8,1,2) = [  0      0       t    1/sqrt(2)   ];
controlPts(:,9,1,2) = [  0      0       t    1           ];

controlPts(:,1,2,2) = [  R      0       t    1           ];
controlPts(:,2,2,2) = [  R      R       t    1/sqrt(2)   ];
controlPts(:,3,2,2) = [  0   	R       t    1           ];
controlPts(:,4,2,2) = [ -R      R       t    1/sqrt(2)   ];
controlPts(:,5,2,2) = [ -R      0       t    1           ];
controlPts(:,6,2,2) = [ -R     -R       t    1/sqrt(2)   ];
controlPts(:,7,2,2) = [  0     -R       t    1           ];
controlPts(:,8,2,2) = [  R     -R       t    1/sqrt(2)   ];
controlPts(:,9,2,2) = [  R  	0       t    1           ];

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});