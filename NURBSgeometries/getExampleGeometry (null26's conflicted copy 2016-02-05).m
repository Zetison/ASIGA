function nurbs = getExampleGeometry()

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 1 1 1];
controlPts = zeros(3,4,3);
R_i = 1;
R_m = 2;
R_o = 3;
tw1 = 0.9;
tw2 = -0.5;
te = 0.4;

controlPts(1:2,1,1) = [R_i-te        0];
controlPts(1:2,2,1) = [R_i        R_i];
controlPts(1:2,3,1) = [0        R_i];
controlPts(1:2,4,1) = [-2*R_i     R_i];
controlPts(1:2,5,1) = [-2*R_i     0];
controlPts(1:2,6,1) = [-2*R_i    -R_i];
controlPts(1:2,7,1) = [0       -R_i];
controlPts(1:2,8,1) = [R_i       -R_i];
controlPts(1:2,9,1) = [R_i-te        0];

controlPts(1:2,1,2) = [R_i+0.1        tw2];
controlPts(1:2,2,2) = [R_m-tw2        R_m];
controlPts(1:2,3,2) = [-tw2        R_m];
controlPts(1:2,4,2) = [-2*R_m-tw2     R_m];
controlPts(1:2,5,2) = [-2*R_m     -tw2];
controlPts(1:2,6,2) = [-2*R_m+tw2    -R_m];
controlPts(1:2,7,2) = [0       -R_m];
controlPts(1:2,8,2) = [R_m       -R_m+tw2];
controlPts(1:2,9,2) = [R_i+0.1        tw2];

controlPts(1:2,1,3) = [R_o        tw1];
controlPts(1:2,2,3) = [R_o-tw1        R_o];
controlPts(1:2,3,3) = [-tw1        R_o];
controlPts(1:2,4,3) = [-2*R_o-tw1     R_o];
controlPts(1:2,5,3) = [-2*R_o     -tw1];
controlPts(1:2,6,3) = [-2*R_o+tw1    -R_o];
controlPts(1:2,7,3) = [+tw1       -R_o];
controlPts(1:2,8,3) = [R_o       -R_o+tw1];
controlPts(1:2,9,3) = [R_o        tw1];


controlPts(3,:,:) = ones(1,size(controlPts,2),size(controlPts,3));
nurbs = createNURBSobject(controlPts,{Xi, Eta});

nurbs = insertKnotsInNURBS(nurbs,{[] [0.5 0.5]});