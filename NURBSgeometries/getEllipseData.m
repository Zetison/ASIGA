function solid = getEllipseData(c_x,c_y)

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;

controlPts = zeros(3,9);

% outer surface
controlPts(:,1) = [ c_x   0     1           ];
controlPts(:,2) = [ c_x   c_y   1/sqrt(2)   ];
controlPts(:,3) = [ 0   c_y     1           ];
controlPts(:,4) = [ -c_x   c_y     1/sqrt(2)   ];
controlPts(:,5) = [ -c_x   0     1           ];
controlPts(:,6) = [ -c_x   -c_y    1/sqrt(2)   ];
controlPts(:,7) = [ 0   -c_y    1           ];
controlPts(:,8) = [ c_x   -c_y    1/sqrt(2)   ];
controlPts(:,9) = [ c_x   0    1           ];

solid = createNURBSobject(controlPts,Xi);
