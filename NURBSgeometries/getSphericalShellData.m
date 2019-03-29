function solid = getSphericalShellData(R_o, t, alignWithAxis)
if nargin < 3
    alignWithAxis = 'Zaxis';
end
R_i = R_o-t;
Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 1 1 2 2 2]/2;
Zeta = [0 0 1 1];

controlPts = zeros(4,9,5,2);

% inner surface
controlPts(:,1,1,1) = [ 0     0      -R_i   1           ];
controlPts(:,2,1,1) = [ 0     0      -R_i   1/sqrt(2)   ];
controlPts(:,3,1,1) = [ 0     0      -R_i   1           ];
controlPts(:,4,1,1) = [ 0     0      -R_i   1/sqrt(2)   ];
controlPts(:,5,1,1) = [ 0     0      -R_i   1           ];
controlPts(:,6,1,1) = [ 0     0      -R_i   1/sqrt(2)   ];
controlPts(:,7,1,1) = [ 0     0      -R_i   1           ];
controlPts(:,8,1,1) = [ 0     0      -R_i   1/sqrt(2)   ];
controlPts(:,9,1,1) = [ 0     0      -R_i   1           ];

controlPts(:,1,2,1) = [  R_i    0    	-R_i    1/sqrt(2)   ];
controlPts(:,2,2,1) = [  R_i	R_i   	-R_i    1/2         ];
controlPts(:,3,2,1) = [  0      R_i 	-R_i    1/sqrt(2)   ];
controlPts(:,4,2,1) = [ -R_i    R_i   	-R_i    1/2         ];
controlPts(:,5,2,1) = [ -R_i    0    	-R_i    1/sqrt(2)   ];
controlPts(:,6,2,1) = [ -R_i   -R_i  	-R_i    1/2         ];
controlPts(:,7,2,1) = [  0     -R_i   	-R_i    1/sqrt(2)   ];
controlPts(:,8,2,1) = [  R_i   -R_i   	-R_i    1/2         ];
controlPts(:,9,2,1) = [  R_i    0   	-R_i    1/sqrt(2)   ];

controlPts(:,1,3,1) = [  R_i   0    	0	1           ];
controlPts(:,2,3,1) = [  R_i   R_i   	0 	1/sqrt(2)   ];
controlPts(:,3,3,1) = [	 0     R_i  	0  	1           ];
controlPts(:,4,3,1) = [	-R_i   R_i  	0 	1/sqrt(2)   ];
controlPts(:,5,3,1) = [	-R_i   0    	0 	1           ];
controlPts(:,6,3,1) = [	-R_i  -R_i   	0	1/sqrt(2)   ];
controlPts(:,7,3,1) = [	 0    -R_i   	0	1           ];
controlPts(:,8,3,1) = [	 R_i  -R_i   	0	1/sqrt(2)   ];
controlPts(:,9,3,1) = [	 R_i   0    	0	1           ];

controlPts(:,1,4,1) = [  R_i 	0   	R_i     1/sqrt(2)   ];
controlPts(:,2,4,1) = [  R_i 	R_i   	R_i     1/2         ];
controlPts(:,3,4,1) = [  0      R_i  	R_i     1/sqrt(2)   ];
controlPts(:,4,4,1) = [ -R_i    R_i  	R_i     1/2         ];
controlPts(:,5,4,1) = [ -R_i    0    	R_i     1/sqrt(2)   ];
controlPts(:,6,4,1) = [ -R_i   -R_i 	R_i     1/2         ];
controlPts(:,7,4,1) = [  0     -R_i  	R_i     1/sqrt(2)   ];
controlPts(:,8,4,1) = [  R_i   -R_i  	R_i     1/2         ];
controlPts(:,9,4,1) = [  R_i    0       R_i     1/sqrt(2)   ];

controlPts(:,1,5,1) = [  0     0      R_i   1           ];
controlPts(:,2,5,1) = [  0     0      R_i   1/sqrt(2)   ];
controlPts(:,3,5,1) = [  0     0      R_i   1           ];
controlPts(:,4,5,1) = [  0     0      R_i   1/sqrt(2)   ];
controlPts(:,5,5,1) = [  0     0      R_i   1           ];
controlPts(:,6,5,1) = [  0     0      R_i   1/sqrt(2)   ];
controlPts(:,7,5,1) = [  0     0      R_i   1           ];
controlPts(:,8,5,1) = [  0     0      R_i   1/sqrt(2)   ];
controlPts(:,9,5,1) = [  0     0      R_i   1           ];

% outer surface
controlPts(:,1,1,2) = [ 0     0      -R_o   1           ];
controlPts(:,2,1,2) = [ 0     0      -R_o   1/sqrt(2)   ];
controlPts(:,3,1,2) = [ 0     0      -R_o   1           ];
controlPts(:,4,1,2) = [ 0     0      -R_o   1/sqrt(2)   ];
controlPts(:,5,1,2) = [ 0     0      -R_o   1           ];
controlPts(:,6,1,2) = [ 0     0      -R_o   1/sqrt(2)   ];
controlPts(:,7,1,2) = [ 0     0      -R_o   1           ];
controlPts(:,8,1,2) = [ 0     0      -R_o   1/sqrt(2)   ];
controlPts(:,9,1,2) = [ 0     0      -R_o   1           ];

controlPts(:,1,2,2) = [  R_o    0    	-R_o    1/sqrt(2)   ];
controlPts(:,2,2,2) = [  R_o	R_o   	-R_o    1/2         ];
controlPts(:,3,2,2) = [  0      R_o 	-R_o    1/sqrt(2)   ];
controlPts(:,4,2,2) = [ -R_o    R_o   	-R_o    1/2         ];
controlPts(:,5,2,2) = [ -R_o    0    	-R_o    1/sqrt(2)   ];
controlPts(:,6,2,2) = [ -R_o   -R_o  	-R_o    1/2         ];
controlPts(:,7,2,2) = [  0     -R_o   	-R_o    1/sqrt(2)   ];
controlPts(:,8,2,2) = [  R_o   -R_o   	-R_o    1/2         ];
controlPts(:,9,2,2) = [  R_o    0   	-R_o    1/sqrt(2)   ];

controlPts(:,1,3,2) = [  R_o   0    	0	1           ];
controlPts(:,2,3,2) = [  R_o   R_o   	0 	1/sqrt(2)   ];
controlPts(:,3,3,2) = [	 0     R_o  	0  	1           ];
controlPts(:,4,3,2) = [	-R_o   R_o  	0 	1/sqrt(2)   ];
controlPts(:,5,3,2) = [	-R_o   0    	0 	1           ];
controlPts(:,6,3,2) = [	-R_o  -R_o   	0	1/sqrt(2)   ];
controlPts(:,7,3,2) = [	 0    -R_o   	0	1           ];
controlPts(:,8,3,2) = [	 R_o  -R_o   	0	1/sqrt(2)   ];
controlPts(:,9,3,2) = [	 R_o   0    	0	1           ];

controlPts(:,1,4,2) = [  R_o 	0   	R_o     1/sqrt(2)   ];
controlPts(:,2,4,2) = [  R_o 	R_o   	R_o     1/2         ];
controlPts(:,3,4,2) = [  0      R_o  	R_o     1/sqrt(2)   ];
controlPts(:,4,4,2) = [ -R_o    R_o  	R_o     1/2         ];
controlPts(:,5,4,2) = [ -R_o    0    	R_o     1/sqrt(2)   ];
controlPts(:,6,4,2) = [ -R_o   -R_o 	R_o     1/2         ];
controlPts(:,7,4,2) = [  0     -R_o  	R_o     1/sqrt(2)   ];
controlPts(:,8,4,2) = [  R_o   -R_o  	R_o     1/2         ];
controlPts(:,9,4,2) = [  R_o    0       R_o     1/sqrt(2)   ];

controlPts(:,1,5,2) = [  0     0      R_o   1           ];
controlPts(:,2,5,2) = [  0     0      R_o   1/sqrt(2)   ];
controlPts(:,3,5,2) = [  0     0      R_o   1           ];
controlPts(:,4,5,2) = [  0     0      R_o   1/sqrt(2)   ];
controlPts(:,5,5,2) = [  0     0      R_o   1           ];
controlPts(:,6,5,2) = [  0     0      R_o   1/sqrt(2)   ];
controlPts(:,7,5,2) = [  0     0      R_o   1           ];
controlPts(:,8,5,2) = [  0     0      R_o   1/sqrt(2)   ];
controlPts(:,9,5,2) = [  0     0      R_o   1           ];

switch alignWithAxis
    case 'Xaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = temp;
    case 'Yaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = temp;
    case 'Zaxis'
        % Nothing to be done
end



solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
