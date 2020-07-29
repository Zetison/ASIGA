function nurbs = getTorusPartData(R,t)
error('Depricated')
Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

controlPts = zeros(4,3,3);
a = sqrt(R^2-(t/2)^2);
phi = acos(t/2/a);
phi2 = 2*asin(t/2/R);
w = cos(phi/2);
w2 = cos(phi2/2);
% b = a-t/a*(R-t);
a2 = sqrt(a^2-(t/2)^2);
b = a2-t/2/a2*(a2-t/2);
c = R/cos(phi2/2);

controlPts(:,1,1) = [ 0   a    -t/2      1           ];
controlPts(:,2,1) = [ 0   c     0        w2            ];
controlPts(:,3,1) = [ 0   a     t/2      1           ];

controlPts(:,1,2) = [ b  a     -t/2        w           ];
controlPts(:,2,2) = [ b  a     0        w2*w       ];
controlPts(:,3,2) = [ b  a     t/2        w           ];

controlPts(:,1,3) = [ a2  t/2     -t/2        1           ];
controlPts(:,2,3) = [ a2  t/2     0        w2       ];
controlPts(:,3,3) = [ a2  t/2     t/2        1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});


