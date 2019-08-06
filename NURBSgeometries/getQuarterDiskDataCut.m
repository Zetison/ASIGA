function nurbs = getQuarterDiskDataCut(R,t)

Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];

a = sqrt(R^2-t^2);

controlPts = zeros(4,3,2);
phi = acos(t/a);
w = cos(phi/2);
a2 = sqrt(a^2-t^2);
b = a2-t/a2*(a2-t);

controlPts(:,1,1) = [ t   0     0      1           ];
controlPts(:,2,1) = [ t   0     0      w            ];
controlPts(:,3,1) = [ t   0     0      1           ];

controlPts(:,1,2) = [ a  0     0        1           ];
controlPts(:,2,2) = [ a  b     0        w       ];
controlPts(:,3,2) = [ t  a2     0        1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});


