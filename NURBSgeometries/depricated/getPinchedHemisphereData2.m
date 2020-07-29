function solid = getPinchedHemisphereData2(R_i)
error('Depricated')


Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];
controlPts = zeros(4,3,3);

controlPts(:,1,1) = [  R_i	0       0       1            ];
controlPts(:,2,1) = [  R_i	R_i     0       1/sqrt(2)	 ];
controlPts(:,3,1) = [  0      R_i     0       1            ];
controlPts(:,1,2) = [  R_i	0       R_i     1/sqrt(2)    ];
controlPts(:,2,2) = [  R_i	R_i     R_i     1/2          ];
controlPts(:,3,2) = [  0      R_i     R_i     1/sqrt(2)    ];
controlPts(:,1,3) = [  0      0       R_i     1            ];
controlPts(:,2,3) = [  0      0   	R_i     1/sqrt(2)    ];
controlPts(:,3,3) = [  0      0   	R_i     1            ];
controlPts = permute(controlPts,[1,3,2]);
solid = createNURBSobject(controlPts,{Xi, Eta});