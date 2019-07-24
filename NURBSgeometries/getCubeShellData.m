function nurbs = getCubeShellData(s,t)

Xi = [0 0 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

nurbs = cell(1,6); 
controlPts = zeros(4,2,2,2); 

controlPts(:,1,1,1) = [0+t, 0+t, s-t, 1];
controlPts(:,2,1,1) = [s-t, 0+t, s-t, 1];
controlPts(:,1,2,1) = [0+t, s-t, s-t, 1];
controlPts(:,2,2,1) = [s-t, s-t, s-t, 1];

controlPts(:,1,1,2) = [0, 0, s, 1];
controlPts(:,2,1,2) = [s, 0, s, 1];
controlPts(:,1,2,2) = [0, s, s, 1];
controlPts(:,2,2,2) = [s, s, s, 1];
 
nurbs{1} = createNURBSobject(controlPts,{Xi, Eta, Zeta}); 
nurbs(1) = translateNURBS(nurbs(1),-[s/2,s/2,s/2]);

for i = 2:4
    nurbs(i) = rotateNURBS(nurbs{1},(i-1)*pi/2,'Yaxis');
end
nurbs(5) = rotateNURBS(nurbs{1},pi/2,'Xaxis');
nurbs(6) = rotateNURBS(nurbs{1},-pi/2,'Xaxis'); 
nurbs = translateNURBS(nurbs,[s/2,s/2,s/2]);