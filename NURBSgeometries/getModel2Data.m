function nurbs = getModel2Data(R,t,L,phi)
% % 
nurbs{1} = getQuarterDiskData2(R);
nurbs(2) = rotateNURBS(nurbs(1),-pi/2,'Yaxis');
nurbs(1:2) = rotateNURBS(nurbs(1:2),pi/2,'Xaxis');
nurbs(3) = rotateNURBS(nurbs(1),pi/2,'Xaxis');
nurbs(1:3) = elevateDegreeInPatches(nurbs(1:3),[0 1]);
nurbs{4} = getSphericalShellDataPatchedQuarter2(R);
nurbs(4) = rotateNURBS(nurbs(4),-pi/2,'Xaxis');
nurbs(5) = rotateNURBS(nurbs(4),pi,'Xaxis');
nurbs(6) = rotateNURBS(nurbs(4),pi/2,'Xaxis');

nurbs(1:6) = translateNURBS(nurbs(1:6),[L/2,0,0]);

nurbs{7} = getCylinderData(R, L);
nurbs(7) = rotateNURBS(nurbs(7),pi/2,'Yaxis');
nurbs(7) = translateNURBS(nurbs(7),[-L/2,0,0]);
nurbs(7) = elevateDegreeInPatches(nurbs(7),[0 1]);

nurbs{8} = getQuarterDiskDataCut(R,t/2);
nurbs(8) = mirrorNURBS(nurbs(8),'xy');
nurbs(8) = translateNURBS(nurbs(8),[0,-t/2,0]);
nurbs(9) = rotateNURBS(nurbs(8),pi/2,'Xaxis');
nurbs{10} = getQuarterDiskDataCut2(R,t/2);
nurbs(10) = rotateNURBS(nurbs(10),-pi/2,'Yaxis');
nurbs(8:10) = elevateDegreeInPatches(nurbs(8:10),[0 1]);
nurbs{11} = getTorusPartData(R,t);
nurbs(8:9) = translateNURBS(nurbs(8:9),[0,t/2,t/2]);
nurbs(12:15) = rotateNURBS(nurbs(8:11),pi/2,'Xaxis');
nurbs(16:19) = rotateNURBS(nurbs(8:11),pi,'Xaxis');
nurbs(20:23) = rotateNURBS(nurbs(8:11),3*pi/2,'Xaxis');


Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

controlPts = zeros(4,3,3);
% phi = pi/2-2*asin(t/R);
% w = cos(phi/2);
a = sqrt(R^2-(t/2)^2);
a2 = sqrt(a^2-(t/2)^2);
% b = (a^2+t^2)/(a+t);
phi2 = 2*asin(t/2/R);
c = R/cos(phi2/2);
w2 = cos(phi2/2);

controlPts(:,:,1) = nurbs{11}.coeffs(:,:,end);
controlPts(:,:,3) = nurbs{19}.coeffs(:,3:-1:1,end);
controlPts(:,1,2) = nurbs{23}.coeffs(:,2,end);
controlPts(:,3,2) = nurbs{15}.coeffs(:,2,end);
controlPts(:,2,2) = [a2,0,0,w2^2];

nurbs{24} = createNURBSobject(controlPts,{Xi, Eta});

nurbs(8:24) = rotateNURBS(nurbs(8:24),pi,'Zaxis');
nurbs(8:24) = translateNURBS(nurbs(8:24),[-L/2,0,0]);
nurbs = rotateNURBS(nurbs,phi,'Xaxis');
