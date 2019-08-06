function nurbs = getQuarterDiskData2(R)

Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
controlPts = zeros(4,3,2);

controlPts(:,1,1) = [ 0   0     0      1           ];
controlPts(:,2,1) = [ 0   0     0      1/sqrt(2)   ];
controlPts(:,3,1) = [ 0   0     0      1           ];

controlPts(:,1,2) = [ R  0     0        1           ];
controlPts(:,2,2) = [ R  R     0        1/sqrt(2)	];
controlPts(:,3,2) = [ 0  R     0        1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});


