function nurbs = getQuarterDiskDataCut2(R,t)

error('Depricated: use getWedgeData')
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];

controlPts = zeros(4,3,2);
phi = pi/2-2*asin(t/R);
w = cos(phi/2);
a = sqrt(R^2-t^2);
b = (a^2+t^2)/(a+t);

controlPts(:,1,1) = [ t   t     0      1           ];
controlPts(:,2,1) = [ t   t     0      w            ];
controlPts(:,3,1) = [ t   t     0      1           ];

controlPts(:,1,2) = [ a  t     0        1           ];
controlPts(:,2,2) = [ b  b     0        w       ];
controlPts(:,3,2) = [ t  a     0        1           ];

nurbs = createNURBSobject(controlPts,{Xi, Eta});


