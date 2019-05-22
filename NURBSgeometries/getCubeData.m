function nurbs = getCubeData(s)

Xi = [0 0 1 1];
Eta = [0 0 1 1];

nurbs = cell(1,6); 
controlPts = zeros(4,2,2); 

controlPts(:,1,1) = [0, 0, s, 1];
controlPts(:,2,1) = [s, 0, s, 1];
controlPts(:,1,2) = [0, s, s, 1];
controlPts(:,2,2) = [s, s, s, 1];

controlPts(1:3,:,:) = controlPts(1:3,:,:) - s/2;
 
nurbs{1} = createNURBSobject(controlPts,{Xi, Eta}); 

for i = 2:4
    nurbs(i) = rotateNURBS(nurbs{1},(i-1)*pi/2,'Yaxis');
end
nurbs(5) = rotateNURBS(nurbs{1},pi/2,'Xaxis');
nurbs(6) = rotateNURBS(nurbs{1},-pi/2,'Xaxis');