function nurbs = getSphericalShellDataPatchedQuarter2(R)
error('Use getEllipsoidData() instead')
Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

controlPts = zeros(4,3,3);

% inner surface
controlPts(:,1,1) = [ 0     0      -R   1           ];
controlPts(:,2,1) = [ 0     0      -R   1/sqrt(2)   ];
controlPts(:,3,1) = [ 0     0      -R   1           ];

controlPts(:,1,2) = [  R    0    	-R    1/sqrt(2)   ];
controlPts(:,2,2) = [  R	R   	-R    1/2         ];
controlPts(:,3,2) = [  0      R 	-R    1/sqrt(2)   ];

controlPts(:,1,3) = [  R   0    	0	1           ];
controlPts(:,2,3) = [  R   R   	0 	1/sqrt(2)   ];
controlPts(:,3,3) = [  0     R  	0  	1           ];
% 
% controlPts(:,1,4) = [  R 	0   	R     1/sqrt(2)   ];
% controlPts(:,2,4) = [  R 	R   	R     1/2         ];
% controlPts(:,3,4) = [  0      R  	R     1/sqrt(2)   ];
% controlPts(:,4,4) = [ -R    R  	R     1/2         ];
% controlPts(:,5,4) = [ -R    0    	R     1/sqrt(2)   ];
% controlPts(:,6,4) = [ -R   -R 	R     1/2         ];
% controlPts(:,7,4) = [  0     -R  	R     1/sqrt(2)   ];
% controlPts(:,8,4) = [  R   -R  	R     1/2         ];
% controlPts(:,9,4) = [  R    0       R     1/sqrt(2)   ];
% 
% controlPts(:,1,5) = [  0     0      R   1           ];
% controlPts(:,2,5) = [  0     0      R   1/sqrt(2)   ];
% controlPts(:,3,5) = [  0     0      R   1           ];
% controlPts(:,4,5) = [  0     0      R   1/sqrt(2)   ];
% controlPts(:,5,5) = [  0     0      R   1           ];
% controlPts(:,6,5) = [  0     0      R   1/sqrt(2)   ];
% controlPts(:,7,5) = [  0     0      R   1           ];
% controlPts(:,8,5) = [  0     0      R   1/sqrt(2)   ];
% controlPts(:,9,5) = [  0     0      R   1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});
