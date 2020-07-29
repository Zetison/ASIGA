function solid = getStarShapeData_2D()
error('Depricated')

Xi = [0 0 1 2 3 4 5 6 7 8 8]/8;

controlPts = zeros(3,9);

% outer surface
controlPts(:,1) = [ 1       0       1           ];
controlPts(:,2) = [ 0.25    0.25    1           ];
controlPts(:,3) = [ 0       1       1           ];
controlPts(:,4) = [ -0.25   0.25    1           ];
controlPts(:,5) = [ -1      0       1           ];
controlPts(:,6) = [ -0.25   -0.25    1           ];
controlPts(:,7) = [ 0       -1       1           ];
controlPts(:,8) = [ 0.25    -0.25    1           ];
controlPts(:,9) = [ 1       0       1           ];

solid = createNURBSobject(controlPts,Xi);
