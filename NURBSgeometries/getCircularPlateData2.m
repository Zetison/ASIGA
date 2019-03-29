function solid = getCircularPlateData2(R, t)

Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,3,3,2);

controlPts(:,1,1,1) = [  R          0           0       1           ];
controlPts(:,2,1,1) = [  R	       -R        	0       1/sqrt(2)   ];
controlPts(:,3,1,1) = [  0         -R           0       1           ];

controlPts(:,1,2,1) = [  R          R           0       1/sqrt(2)   ];
controlPts(:,2,2,1) = [  0          0        	0       sqrt(2)-1   ];
controlPts(:,3,2,1) = [ -R         -R           0       1/sqrt(2)   ];

controlPts(:,1,3,1) = [  0          R           0       1           ];
controlPts(:,2,3,1) = [ -R          R        	0       1/sqrt(2)   ];
controlPts(:,3,3,1) = [ -R          0           0       1           ];

controlPts(:,1,1,2) = [  R          0           t       1           ];
controlPts(:,2,1,2) = [  R         -R        	t       1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0         -R           t       1           ];

controlPts(:,1,2,2) = [  R          R           t       1/sqrt(2)   ];
controlPts(:,2,2,2) = [  0          0        	t       sqrt(2)-1   ];
controlPts(:,3,2,2) = [ -R         -R           t       1/sqrt(2)   ];

controlPts(:,1,3,2) = [  0          R           t       1           ];
controlPts(:,2,3,2) = [ -R          R        	t       1/sqrt(2)   ];
controlPts(:,3,3,2) = [ -R          0           t       1           ];

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});