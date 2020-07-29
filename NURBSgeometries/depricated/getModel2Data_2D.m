function solid = getModel2Data_2D(R1, R2, L, t, x_0)

error('Depricated use getModelBeTSSiM2Data')

xi1 = 0.14;
xi2 = 0.37;
xi3 = 0.4;
xi4 = 0.45;
xi5 = 0.55;
xi6 = 0.6;
xi7 = 0.63;
xi8 = 0.86;

Xi = [0 0 0 xi1 xi1 xi2 xi2 xi3 xi3 xi4 xi4 xi5 xi5 xi6 xi6 xi7 xi7 xi8 xi8 1 1 1];

controlPts = zeros(3,9);

mid = (R2+t)/2;

% outer surface
controlPts(:,1) = [ R1       0       1           ];
controlPts(:,2) = [ R1       R1       1/sqrt(2)   ];
controlPts(:,3) = [ 0       R1       1           ];
controlPts(:,4) = [ -L/2    (R1+R2)/2       1           ];
controlPts(:,5) = [ -L      R2       1           ];
controlPts(:,6) = [ -L      mid     1           ];
controlPts(:,7) = [ -L      t       1           ];
controlPts(:,8) = [ -L-R2/2  t       1           ];
controlPts(:,9) = [ -L-R2    t       1           ];
controlPts(:,10) = [ -L-R2    0       1          ];
controlPts(:,11) = [ -L-R2    -t      1          ];
controlPts(:,12) = [ -L-R2/2  -t      1          ];
controlPts(:,13) = [ -L      -t       1         ];
controlPts(:,14) = [ -L      -mid     1         ];
controlPts(:,15) = [ -L      -R2       1         ];
controlPts(:,16) = [ -L/2    -(R1+R2)/2       1         ];
controlPts(:,17) = [ 0       -R1       1         ];
controlPts(:,18) = [ R1       -R1       1/sqrt(2) ];
controlPts(:,19) = [ R1       0        1         ];

controlPts(1,:) = controlPts(1,:) + x_0(1);
controlPts(2,:) = controlPts(2,:) + x_0(2);

solid = createNURBSobject(controlPts,Xi);
