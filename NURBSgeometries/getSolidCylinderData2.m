function solid = getSolidCylinderData2(R_i, R_o, L)

% Quaarter of the full cylinder

Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,3,2,2);

controlPts(:,1,1,1) = [  R_i	0       0	 1           ];
controlPts(:,2,1,1) = [  R_i	R_i     0	 1/sqrt(2)   ];
controlPts(:,3,1,1) = [  0      R_i     0    1           ];

controlPts(:,1,2,1) = [  R_i	0       L    1           ];
controlPts(:,2,2,1) = [  R_i	R_i     L    1/sqrt(2)   ];
controlPts(:,3,2,1) = [  0  	R_i     L    1           ];

controlPts(:,1,1,2) = [  R_o	0       0    1           ];
controlPts(:,2,1,2) = [  R_o	R_o     0    1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0   	R_o     0    1           ];

controlPts(:,1,2,2) = [  R_o	0       L    1           ];
controlPts(:,2,2,2) = [  R_o	R_o     L    1/sqrt(2)   ];
controlPts(:,3,2,2) = [  0   	R_o     L    1           ];

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});