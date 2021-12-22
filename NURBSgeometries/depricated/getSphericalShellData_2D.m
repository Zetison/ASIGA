function solid = getSphericalShellData_2D(R, R_o,alignWithAxis)
error('Use getEllipsoidData() instead')
if nargin < 3
    alignWithAxis = 'Xaxis';
end
Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 1 1];

controlPts = zeros(3,9,2);

% inner surface
controlPts(:,1,1) = [  R    0    	1   ];
controlPts(:,2,1) = [  R	R   	1/sqrt(2)         ];
controlPts(:,3,1) = [  0      R 	1   ];
controlPts(:,4,1) = [ -R    R   1/sqrt(2)         ];
controlPts(:,5,1) = [ -R    0    	1   ];
controlPts(:,6,1) = [ -R   -R  	1/sqrt(2)         ];
controlPts(:,7,1) = [  0     -R   1   ];
controlPts(:,8,1) = [  R   -R   1/sqrt(2)         ];
controlPts(:,9,1) = [  R    0   	1   ];

% outer surface
controlPts(:,1,2) = [  R_o    0    	1   ];
controlPts(:,2,2) = [  R_o	R_o   	1/sqrt(2)         ];
controlPts(:,3,2) = [  0      R_o 	1   ];
controlPts(:,4,2) = [ -R_o    R_o   1/sqrt(2)         ];
controlPts(:,5,2) = [ -R_o    0    	1   ];
controlPts(:,6,2) = [ -R_o   -R_o  	1/sqrt(2)         ];
controlPts(:,7,2) = [  0     -R_o   1   ];
controlPts(:,8,2) = [  R_o   -R_o   1/sqrt(2)         ];
controlPts(:,9,2) = [  R_o    0   	1   ];

switch alignWithAxis
    case 'Xaxis' 
        % Nothing to be done
    case 'Yaxis' 
        temp = controlPts(1,:,:);
        controlPts(1,:,:) = controlPts(2,:,:);
        controlPts(2,:,:) = temp;
end


solid = createNURBSobject(controlPts,{Xi, Eta});
