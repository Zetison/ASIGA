function solid = getEllipsoidalShellData(c_x,c_y,c_z,t,alignWithAxis)
if nargin < 5
    alignWithAxis = 'Zaxis';
end
Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 1 1 2 2 2]/2;
Zeta = [0 0 1 1];
c_xi = c_x-t;
c_yi = c_y-t;
c_zi = c_z-t;

controlPts = zeros(4,9,5,2);

controlPts(:,1,1,1) = [ 0     0      -c_zi   1           ];
controlPts(:,2,1,1) = [ 0     0      -c_zi   1/sqrt(2)   ];
controlPts(:,3,1,1) = [ 0     0      -c_zi   1           ];
controlPts(:,4,1,1) = [ 0     0      -c_zi   1/sqrt(2)   ];
controlPts(:,5,1,1) = [ 0     0      -c_zi   1           ];
controlPts(:,6,1,1) = [ 0     0      -c_zi   1/sqrt(2)   ];
controlPts(:,7,1,1) = [ 0     0      -c_zi   1           ];
controlPts(:,8,1,1) = [ 0     0      -c_zi   1/sqrt(2)   ];
controlPts(:,9,1,1) = [ 0     0      -c_zi   1           ];

controlPts(:,1,2,1) = [  c_xi    0    	-c_zi    1/sqrt(2)   ];
controlPts(:,2,2,1) = [  c_xi	c_yi   	-c_zi    1/2         ];
controlPts(:,3,2,1) = [  0      c_yi 	-c_zi    1/sqrt(2)   ];
controlPts(:,4,2,1) = [ -c_xi    c_yi   -c_zi    1/2         ];
controlPts(:,5,2,1) = [ -c_xi    0    	-c_zi    1/sqrt(2)   ];
controlPts(:,6,2,1) = [ -c_xi   -c_yi  	-c_zi    1/2         ];
controlPts(:,7,2,1) = [  0     -c_yi   	-c_zi    1/sqrt(2)   ];
controlPts(:,8,2,1) = [  c_xi   -c_yi   -c_zi    1/2         ];
controlPts(:,9,2,1) = [  c_xi    0   	-c_zi    1/sqrt(2)   ];

controlPts(:,1,3,1) = [  c_xi   0    	0	1           ];
controlPts(:,2,3,1) = [  c_xi   c_yi   	0 	1/sqrt(2)   ];
controlPts(:,3,3,1) = [	 0     c_yi  	0  	1           ];
controlPts(:,4,3,1) = [	-c_xi   c_yi  	0 	1/sqrt(2)   ];
controlPts(:,5,3,1) = [	-c_xi   0    	0 	1           ];
controlPts(:,6,3,1) = [	-c_xi  -c_yi   	0	1/sqrt(2)   ];
controlPts(:,7,3,1) = [	 0    -c_yi   	0	1           ];
controlPts(:,8,3,1) = [	 c_xi  -c_yi   	0	1/sqrt(2)   ];
controlPts(:,9,3,1) = [	 c_xi   0    	0	1           ];

controlPts(:,1,4,1) = [  c_xi 	0   	c_zi     1/sqrt(2)   ];
controlPts(:,2,4,1) = [  c_xi 	c_yi   	c_zi     1/2         ];
controlPts(:,3,4,1) = [  0      c_yi  	c_zi     1/sqrt(2)   ];
controlPts(:,4,4,1) = [ -c_xi    c_yi  	c_zi     1/2         ];
controlPts(:,5,4,1) = [ -c_xi    0    	c_zi     1/sqrt(2)   ];
controlPts(:,6,4,1) = [ -c_xi   -c_yi 	c_zi     1/2         ];
controlPts(:,7,4,1) = [  0     -c_yi  	c_zi     1/sqrt(2)   ];
controlPts(:,8,4,1) = [  c_xi   -c_yi  	c_zi     1/2         ];
controlPts(:,9,4,1) = [  c_xi    0       c_zi    1/sqrt(2)   ];

controlPts(:,1,5,1) = [  0     0      c_zi   1           ];
controlPts(:,2,5,1) = [  0     0      c_zi   1/sqrt(2)   ];
controlPts(:,3,5,1) = [  0     0      c_zi   1           ];
controlPts(:,4,5,1) = [  0     0      c_zi   1/sqrt(2)   ];
controlPts(:,5,5,1) = [  0     0      c_zi   1           ];
controlPts(:,6,5,1) = [  0     0      c_zi   1/sqrt(2)   ];
controlPts(:,7,5,1) = [  0     0      c_zi   1           ];
controlPts(:,8,5,1) = [  0     0      c_zi   1/sqrt(2)   ];
controlPts(:,9,5,1) = [  0     0      c_zi   1           ];

% outer surface
controlPts(:,1,1,2) = [ 0     0      -c_z   1           ];
controlPts(:,2,1,2) = [ 0     0      -c_z   1/sqrt(2)   ];
controlPts(:,3,1,2) = [ 0     0      -c_z   1           ];
controlPts(:,4,1,2) = [ 0     0      -c_z   1/sqrt(2)   ];
controlPts(:,5,1,2) = [ 0     0      -c_z   1           ];
controlPts(:,6,1,2) = [ 0     0      -c_z   1/sqrt(2)   ];
controlPts(:,7,1,2) = [ 0     0      -c_z   1           ];
controlPts(:,8,1,2) = [ 0     0      -c_z   1/sqrt(2)   ];
controlPts(:,9,1,2) = [ 0     0      -c_z   1           ];

controlPts(:,1,2,2) = [  c_x    0    	-c_z    1/sqrt(2)   ];
controlPts(:,2,2,2) = [  c_x	c_y   	-c_z    1/2         ];
controlPts(:,3,2,2) = [  0      c_y 	-c_z    1/sqrt(2)   ];
controlPts(:,4,2,2) = [ -c_x    c_y   	-c_z    1/2         ];
controlPts(:,5,2,2) = [ -c_x    0    	-c_z    1/sqrt(2)   ];
controlPts(:,6,2,2) = [ -c_x   -c_y  	-c_z    1/2         ];
controlPts(:,7,2,2) = [  0     -c_y   	-c_z    1/sqrt(2)   ];
controlPts(:,8,2,2) = [  c_x   -c_y   	-c_z    1/2         ];
controlPts(:,9,2,2) = [  c_x    0   	-c_z    1/sqrt(2)   ];

controlPts(:,1,3,2) = [  c_x   0    	0	1           ];
controlPts(:,2,3,2) = [  c_x   c_y   	0 	1/sqrt(2)   ];
controlPts(:,3,3,2) = [	 0     c_y  	0  	1           ];
controlPts(:,4,3,2) = [	-c_x   c_y  	0 	1/sqrt(2)   ];
controlPts(:,5,3,2) = [	-c_x   0    	0 	1           ];
controlPts(:,6,3,2) = [	-c_x  -c_y   	0	1/sqrt(2)   ];
controlPts(:,7,3,2) = [	 0    -c_y   	0	1           ];
controlPts(:,8,3,2) = [	 c_x  -c_y   	0	1/sqrt(2)   ];
controlPts(:,9,3,2) = [	 c_x   0    	0	1           ];

controlPts(:,1,4,2) = [  c_x 	0   	c_z     1/sqrt(2)   ];
controlPts(:,2,4,2) = [  c_x 	c_y   	c_z     1/2         ];
controlPts(:,3,4,2) = [  0      c_y  	c_z     1/sqrt(2)   ];
controlPts(:,4,4,2) = [ -c_x    c_y  	c_z     1/2         ];
controlPts(:,5,4,2) = [ -c_x    0    	c_z     1/sqrt(2)   ];
controlPts(:,6,4,2) = [ -c_x   -c_y 	c_z     1/2         ];
controlPts(:,7,4,2) = [  0     -c_y  	c_z     1/sqrt(2)   ];
controlPts(:,8,4,2) = [  c_x   -c_y  	c_z     1/2         ];
controlPts(:,9,4,2) = [  c_x    0       c_z     1/sqrt(2)   ];

controlPts(:,1,5,2) = [  0     0      c_z   1           ];
controlPts(:,2,5,2) = [  0     0      c_z   1/sqrt(2)   ];
controlPts(:,3,5,2) = [  0     0      c_z   1           ];
controlPts(:,4,5,2) = [  0     0      c_z   1/sqrt(2)   ];
controlPts(:,5,5,2) = [  0     0      c_z   1           ];
controlPts(:,6,5,2) = [  0     0      c_z   1/sqrt(2)   ];
controlPts(:,7,5,2) = [  0     0      c_z   1           ];
controlPts(:,8,5,2) = [  0     0      c_z   1/sqrt(2)   ];
controlPts(:,9,5,2) = [  0     0      c_z   1           ];

% inner surface
switch alignWithAxis
    case 'Xaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = c_x/c_z*controlPts(3,:,:,:);
        controlPts(3,:,:,:) = c_z/c_y*controlPts(2,:,:,:);
        controlPts(2,:,:,:) = c_y/c_x*temp;
    case 'Yaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = c_x/c_y*controlPts(2,:,:,:);
        controlPts(2,:,:,:) = c_y/c_z*controlPts(3,:,:,:);
        controlPts(3,:,:,:) = c_z/c_x*temp;
    case 'Zaxis'
        % Nothing to be done
end


solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
