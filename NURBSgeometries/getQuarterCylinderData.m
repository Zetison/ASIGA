function nurbs = getQuarterCylinderData(R,H)

Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];
controlPts = zeros(4,3,2,2);

controlPts(:,1,1,1) = [ 0   0     0      1           ];
controlPts(:,2,1,1) = [ 0   0     0      1/sqrt(2)   ];
controlPts(:,3,1,1) = [ 0   0     0      1           ];

controlPts(:,1,2,1) = [ R  0     0        1           ];
controlPts(:,2,2,1) = [ R  R     0        1/sqrt(2)	];
controlPts(:,3,2,1) = [ 0  R     0        1           ];

controlPts(:,1,1,2) = [ 0   0     H      1           ];
controlPts(:,2,1,2) = [ 0   0     H      1/sqrt(2)   ];
controlPts(:,3,1,2) = [ 0   0     H      1           ];

controlPts(:,1,2,2) = [ R  0     H        1           ];
controlPts(:,2,2,2) = [ R  R     H        1/sqrt(2)	];
controlPts(:,3,2,2) = [ 0  R     H        1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta, Zeta});