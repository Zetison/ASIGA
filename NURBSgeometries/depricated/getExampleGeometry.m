function nurbs = getExampleGeometry()
error('Depricated')
Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 0.5 0.5 1 1 1];
controlPts = zeros(3,4,5);

c_x = 5;
c_y = 0.8;
c_xi = 0.9*c_x;
c_yi = 0.9*c_y;
s = 0.9;
controlPts(3, 1:2:end, [1, 3]) = 1; 
controlPts(3,2:2:end,[1, 3]) = 1/sqrt(2);


controlPts(1:3,1,3) = [c_xi      0       1];
controlPts(1:3,2,3) = [c_xi     c_yi     1/sqrt(2)];
controlPts(1:3,3,3) = [0        c_yi/10     1];
controlPts(1:3,4,3) = [-c_xi 	c_yi     1/sqrt(2)];
controlPts(1:3,5,3) = [-c_xi     0       1];
controlPts(1:3,6,3) = [-c_xi    -c_yi/5     1/sqrt(2)];
controlPts(1:3,7,3) = [0       -c_yi/2     1];
controlPts(1:3,8,3) = [c_xi     -c_yi     1/sqrt(2)];
controlPts(1:3,9,3) = [c_xi      0       1];

controlPts(1:3,1,5) = [c_x      0       1];
controlPts(1:3,2,5) = [c_x      c_y     1/sqrt(2)];
controlPts(1:3,3,5) = [0        c_y     1];
controlPts(1:3,4,5) = [-c_x 	c_y     1/sqrt(2)];
controlPts(1:3,5,5) = [-c_x     0       1];
controlPts(1:3,6,5) = [-c_x    -c_y     1/sqrt(2)];
controlPts(1:3,7,5) = [0       -c_y     1];
controlPts(1:3,8,5) = [c_x     -c_y     1/sqrt(2)];
controlPts(1:3,9,5) = [c_x      0       1];

for i = 1:9
    controlPts(1:2,i,2) = 0.5*controlPts(1:2,i,1) + 0.5*controlPts(1:2,i,3);
end

for i = 1:9
    controlPts(1:2,i,4) = s*controlPts(1:2,i,3) + (1-s)*controlPts(1:2,i,5);
end


controlPts(3,:,:) = ones(1,size(controlPts,2),size(controlPts,3));
nurbs = createNURBSobject(controlPts,{Xi, Eta});

nurbs = insertKnotsInNURBS(nurbs,{[] [0.5 0.5]});