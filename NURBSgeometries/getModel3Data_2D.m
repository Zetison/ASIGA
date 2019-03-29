function solid = getModel3Data_2D(R1, R2, L, x_0)

xi1 = 0.1;
xi2 = 0.37;
xi4 = 1-xi2;
xi5 = 1-xi1;

Xi = [0 0 0 xi1 xi1 xi2 xi2 0.5 0.5 xi4 xi4 xi5 xi5  1 1 1];

controlPts = zeros(3,13);

% outer surface
controlPts(:,1) = [ R1       0       1           ];
controlPts(:,2) = [ R1       R1       1/sqrt(2)   ];
controlPts(:,3) = [ 0       R1       1           ];
controlPts(:,4) = [ -L/2    (R1+R2)/2       1           ];
controlPts(:,5) = [ -L      R2       1           ];

controlPts(:,6) = [ -L-R2      R2     1/sqrt(2)           ];
controlPts(:,7) = [ -L-R2      0       1           ];
controlPts(:,8) = [ -L-R2    -R2       1/sqrt(2)           ];

controlPts(:,9) = [ -L      -R2       1         ];
controlPts(:,10) = [ -L/2    -(R1+R2)/2       1         ];
controlPts(:,11) = [ 0       -R1       1         ];
controlPts(:,12) = [ R1       -R1       1/sqrt(2) ];
controlPts(:,13) = [ R1       0        1         ];

controlPts(1,:) = controlPts(1,:) + x_0(1);
controlPts(2,:) = controlPts(2,:) + x_0(2);

solid = createNURBSobject(controlPts,Xi);
