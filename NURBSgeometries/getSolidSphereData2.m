function solid = getSolidSphereData2(R,alignWithAxis)%,www,p)
if nargin < 2
    alignWithAxis = 'Zaxis';
end
Xi = [0 0 0 0 0 1 1 1 1 1];
Eta = [0 0 0 0 0 1 1 1 1 1];
Zeta = [0 0 0 0 0 1 1 1 1 1];


sr2 = sqrt(2);                                                        
sr3 = sqrt(3);                                                        
sr6 = sqrt(6);   
controlPts = zeros(4,5,5,5);    
controlPts(4,:,:,:) = 1;                 
controlPts(:,1,1,1) = [-4*(sr3-1)   4*(1-sr3)       4*(1-sr3)       4*(3-sr3)];
controlPts(:,2,1,1) = [-sr2         sr2*(sr3-4)     sr2*(sr3-4)     sr2*(3*sr3-2)];
controlPts(:,3,1,1) = [0            4/3*(1-2*sr3)   4/3*(1-2*sr3)	4/3*(5-sr3)];
controlPts(:,4,1,1) = [sr2          sr2*(sr3-4)     sr2*(sr3-4)     sr2*(3*sr3-2)];
controlPts(:,5,1,1) = [4*(sr3-1)	4*(1-sr3)       4*(1-sr3)       4*(3-sr3)];

controlPts(:,1,2,1) = [-sr2*(4-sr3)	-sr2            sr2*(sr3-4) sr2*(3*sr3-2)];
controlPts(:,2,2,1) = [-(3*sr3-2)/2	(2-3*sr3)/2     -(sr3+6)/2  (sr3+6)/2];
controlPts(:,3,2,1) = [0            sr2*(2*sr3-7)/3	-5*sr6/3    sr2*(sr3+6)/3];
controlPts(:,4,2,1) = [(3*sr3-2)/2	(2-3*sr3)/2     -(sr3+6)/2  (sr3+6)/2];
controlPts(:,5,2,1) = [sr2*(4-sr3)	-sr2            sr2*(sr3-4) sr2*(3*sr3-2)];

controlPts(:,1,3,1) = [-4/3*(2*sr3-1)	0           4/3*(1-2*sr3)   4*(5-sr3)/3];
controlPts(:,2,3,1) = [-sr2/3*(7-2*sr3)	0           -5*sr6/3        sr2*(sr3+6)/3];
controlPts(:,3,3,1) = [0                0           4*(sr3-5)/3     4*(5*sr3-1)/9];
controlPts(:,4,3,1) = [sr2/3*(7-2*sr3)	0           -5*sr6/3        sr2*(sr3+6)/3];
controlPts(:,5,3,1) = [4/3*(2*sr3-1)  	0           4/3*(1-2*sr3)   4*(5-sr3)/3];

controlPts(:,1,4,1) = [-sr2*(4-sr3) 	sr2                 sr2*(sr3-4)	sr2*(3*sr3-2)];
controlPts(:,2,4,1) = [-(3*sr3-2)/2     -(2-3*sr3)/2        -(sr3+6)/2	(sr3+6)/2];
controlPts(:,3,4,1) = [0                -sr2*(2*sr3-7)/3 	-5*sr6/3    sr2*(sr3+6)/3];
controlPts(:,4,4,1) = [(3*sr3-2)/2      -(2-3*sr3)/2        -(sr3+6)/2  (sr3+6)/2];
controlPts(:,5,4,1) = [sr2*(4-sr3)  	sr2                 sr2*(sr3-4) sr2*(3*sr3-2)];

controlPts(:,1,5,1) = [-4*(sr3-1) 	-4*(1-sr3)      4*(1-sr3)       4*(3-sr3)];
controlPts(:,2,5,1) = [-sr2         -sr2*(sr3-4)    sr2*(sr3-4)     sr2*(3*sr3-2)];
controlPts(:,3,5,1) = [0            -4/3*(1-2*sr3)  4/3*(1-2*sr3)	4/3*(5-sr3)];
controlPts(:,4,5,1) = [sr2          -sr2*(sr3-4)    sr2*(sr3-4)     sr2*(3*sr3-2)];
controlPts(:,5,5,1) = [4*(sr3-1)  	-4*(1-sr3)      4*(1-sr3)       4*(3-sr3)];

for i = 1:5
    for j = 1:5
        controlPts(1:3,i,j,1) = controlPts(1:3,i,j,1)/controlPts(4,i,j,1);
    end
end  
controlPts(1:3,:,:,1) = controlPts(1:3,:,:,1)*R;  

controlPts(:,:,:,5) = controlPts(:,:,:,1);
controlPts(3,:,:,5) = -controlPts(3,:,:,1);

for i = 2:4
    for j = 1:5
        controlPts(1:3,1,j,6-i) = rotationYaxis(pi/2)*controlPts(1:3,i,j,1);
        controlPts(4,1,j,6-i) = controlPts(4,i,j,1);
        controlPts(1:3,5,j,i) = rotationYaxis(-pi/2)*controlPts(1:3,i,j,1);
        controlPts(4,5,j,i) = controlPts(4,i,j,1);
    end
end   

for i = 2:4
    for j = 2:4
        controlPts(1:3,i,1,6-j) = rotationXaxis(-pi/2)*controlPts(1:3,i,j,1);
        controlPts(4,i,1,6-j) = controlPts(4,i,j,1);
        controlPts(1:3,i,5,j) = rotationXaxis(pi/2)*controlPts(1:3,i,j,1);
        controlPts(4,i,5,j) = controlPts(4,i,j,1);
    end
end   



% Create internal points
x1 = 0.5;
x2 = x1;
x3 = controlPts(3,3,3,5)/2;
w1 = 1;
w2 = 1;
w3 = 1;
w4 = 1;
controlPts(:,2,2,2) = [-x1 -x1 -x1 w1];
controlPts(:,3,2,2) = [0   -x2 -x2 w2];
controlPts(:,4,2,2) = [x1  -x1 -x1 w1];

controlPts(:,2,3,2) = [-x2 0 -x2 w2];
controlPts(:,3,3,2) = [0   0  -x3 w2];
controlPts(:,4,3,2) = [x2 0 -x2 w2];

controlPts(:,2,4,2) = [-x1 x1 -x1 w1];
controlPts(:,3,4,2) = [0   x2 -x2 w2];
controlPts(:,4,4,2) = [x1  x1 -x1 w1];


controlPts(:,2,2,3) = [-x2 -x2 0 w2];
controlPts(:,3,2,3) = [0   -x3 0 w2];
controlPts(:,4,2,3) = [x2  -x2 0 w2];

controlPts(:,2,3,3) = [-x3 0 0 w2];
controlPts(:,3,3,3) = [ 0  0 0 w4];
controlPts(:,4,3,3) = [ x3 0 0 w2];

controlPts(:,2,4,3) = [-x2 x2 0 w2];
controlPts(:,3,4,3) = [0   x3 0 w3];
controlPts(:,4,4,3) = [x2  x2 0 w2];



controlPts(:,2,2,4) = [-x1 -x1 x1 w1];
controlPts(:,3,2,4) = [0   -x2 x2 w2];
controlPts(:,4,2,4) = [x1  -x1 x1 w1];

controlPts(:,2,3,4) = [-x2 0 x2 w2];
controlPts(:,3,3,4) = [0  0  x3 w2];
controlPts(:,4,3,4) = [x2 0 x2 w2];

controlPts(:,2,4,4) = [-x1 x1 x1 w1];
controlPts(:,3,4,4) = [0   x2 x2 w2];
controlPts(:,4,4,4) = [x1  x1 x1 w1];


switch alignWithAxis
    case 'Xaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = temp;
    case 'Yaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = temp;
    case 'Zaxis'
        % Nothing to be done
end



solid = createNURBSobject(controlPts,{Xi Eta Zeta});


function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];


function R_x = rotationYaxis(alpha)

R_x = [cos(alpha),  0,           sin(alpha);
       0,           1,           0;
       -sin(alpha), 0,           cos(alpha)];
