function solid = getPinchedHemisphereData2(R)
error('Depricated')


Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];
controlPts = zeros(4,3,3);

controlPts(:,1,1) = [  R	0       0       1            ];
controlPts(:,2,1) = [  R	R     0       1/sqrt(2)	 ];
controlPts(:,3,1) = [  0      R     0       1            ];
controlPts(:,1,2) = [  R	0       R     1/sqrt(2)    ];
controlPts(:,2,2) = [  R	R     R     1/2          ];
controlPts(:,3,2) = [  0      R     R     1/sqrt(2)    ];
controlPts(:,1,3) = [  0      0       R     1            ];
controlPts(:,2,3) = [  0      0   	R     1/sqrt(2)    ];
controlPts(:,3,3) = [  0      0   	R     1            ];
controlPts = permute(controlPts,[1,3,2]);
solid = createNURBSobject(controlPts,{Xi, Eta});