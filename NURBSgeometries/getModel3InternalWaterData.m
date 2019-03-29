function nurbs = getModel3InternalWaterData(R_i1,R_i2,L,alignWithAxis, eta1, eta2, x_0)
if nargin < 4
    alignWithAxis = 'Zaxis';
end
if nargin < 5
    eta1 = 0.26;
    eta2 = 0.80;
end

if nargin < 7
    x_0 = [0; 0; 0];
end

R_im = (R_i1+R_i2)/2;

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 eta1 eta1 eta2 eta2 1 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,9,7,2);

% inner surface
controlPts(:,1,1,1) = [ -L   0     0      1           ];
controlPts(:,2,1,1) = [ -L   0     0      1/sqrt(2)   ];
controlPts(:,3,1,1) = [ -L   0     0      1           ];
controlPts(:,4,1,1) = [ -L   0     0      1/sqrt(2)   ];
controlPts(:,5,1,1) = [ -L   0     0      1           ];
controlPts(:,6,1,1) = [ -L   0     0      1/sqrt(2)   ];
controlPts(:,7,1,1) = [ -L   0     0      1           ];
controlPts(:,8,1,1) = [ -L   0     0      1/sqrt(2)   ];
controlPts(:,9,1,1) = [ -L   0     0      1           ];

controlPts(:,1,2,1) = [ -L   0	 0    	1/sqrt(2)   ];
controlPts(:,2,2,1) = [ -L   0	 0   	1/2         ];
controlPts(:,3,2,1) = [ -L   0     0 	1/sqrt(2)   ];
controlPts(:,4,2,1) = [ -L   0   0   	1/2         ];
controlPts(:,5,2,1) = [ -L   0   0    	1/sqrt(2)   ];
controlPts(:,6,2,1) = [ -L   0   0  	1/2         ];
controlPts(:,7,2,1) = [ -L   0    0   	1/sqrt(2)   ];
controlPts(:,8,2,1) = [ -L   0   0   	1/2         ];
controlPts(:,9,2,1) = [ -L   0   0   	1/sqrt(2)   ];

controlPts(:,1,3,1) = [ -L     0   0    	1           ];
controlPts(:,2,3,1) = [ -L     0   0   	1/sqrt(2)   ];
controlPts(:,3,3,1) = [	-L     0   0  	1           ];
controlPts(:,4,3,1) = [	-L     0   0  	1/sqrt(2)   ];
controlPts(:,5,3,1) = [	-L     0   0    	1           ];
controlPts(:,6,3,1) = [	-L     0   0   	1/sqrt(2)   ];
controlPts(:,7,3,1) = [	-L     0   0   	1           ];
controlPts(:,8,3,1) = [	-L     0   0   	1/sqrt(2)   ];
controlPts(:,9,3,1) = [	-L     0   0    	1           ];

controlPts(:,1,4,1) = [ -L/2   0   0    	1           ];
controlPts(:,2,4,1) = [ -L/2   0   0   	1/sqrt(2)   ];
controlPts(:,3,4,1) = [ -L/2   0     0  	1           ];
controlPts(:,4,4,1) = [ -L/2   0   0  	1/sqrt(2)   ];
controlPts(:,5,4,1) = [ -L/2   0   0    	1           ];
controlPts(:,6,4,1) = [ -L/2   0   0   	1/sqrt(2)   ];
controlPts(:,7,4,1) = [ -L/2   0     0   	1           ];
controlPts(:,8,4,1) = [ -L/2   0   0   	1/sqrt(2)   ];
controlPts(:,9,4,1) = [ -L/2   0   0    	1           ];

controlPts(:,1,5,1) = [  0     0   0    	1           ];
controlPts(:,2,5,1) = [  0     0   0   	1/sqrt(2)   ];
controlPts(:,3,5,1) = [	 0     0     0  	1           ];
controlPts(:,4,5,1) = [	 0     0   0  	1/sqrt(2)   ];
controlPts(:,5,5,1) = [	 0     0   0    	1           ];
controlPts(:,6,5,1) = [	 0     0   0   	1/sqrt(2)   ];
controlPts(:,7,5,1) = [	 0     0     0   	1           ];
controlPts(:,8,5,1) = [	 0     0   0   	1/sqrt(2)   ];
controlPts(:,9,5,1) = [	 0     0   0    	1           ];

controlPts(:,1,6,1) = [  0   0 	 0   	1/sqrt(2)   ];
controlPts(:,2,6,1) = [  0   0 	 0   	1/2         ];
controlPts(:,3,6,1) = [  0   0     0  	1/sqrt(2)   ];
controlPts(:,4,6,1) = [  0   0 	 0  	1/2         ];
controlPts(:,5,6,1) = [  0   0 	 0    	1/sqrt(2)   ];
controlPts(:,6,6,1) = [  0   0   0      1/2         ];
controlPts(:,7,6,1) = [  0   0     0  	1/sqrt(2)   ];
controlPts(:,8,6,1) = [  0   0   0  	1/2         ];
controlPts(:,9,6,1) = [  0   0   0  	1/sqrt(2)   ];

controlPts(:,1,7,1) = [  0   0     0      1           ];
controlPts(:,2,7,1) = [  0   0     0      1/sqrt(2)   ];
controlPts(:,3,7,1) = [  0   0     0      1           ];
controlPts(:,4,7,1) = [  0   0     0      1/sqrt(2)   ];
controlPts(:,5,7,1) = [  0   0     0      1           ];
controlPts(:,6,7,1) = [  0   0     0      1/sqrt(2)   ];
controlPts(:,7,7,1) = [  0   0     0      1           ];
controlPts(:,8,7,1) = [  0   0     0      1/sqrt(2)   ];
controlPts(:,9,7,1) = [  0   0     0      1           ];

% outer surface
controlPts(:,1,1,2) = [ -L-R_i1   0     0      1           ];
controlPts(:,2,1,2) = [ -L-R_i1   0     0      1/sqrt(2)   ];
controlPts(:,3,1,2) = [ -L-R_i1   0     0      1           ];
controlPts(:,4,1,2) = [ -L-R_i1   0     0      1/sqrt(2)   ];
controlPts(:,5,1,2) = [ -L-R_i1   0     0      1           ];
controlPts(:,6,1,2) = [ -L-R_i1   0     0      1/sqrt(2)   ];
controlPts(:,7,1,2) = [ -L-R_i1   0     0      1           ];
controlPts(:,8,1,2) = [ -L-R_i1   0     0      1/sqrt(2)   ];
controlPts(:,9,1,2) = [ -L-R_i1   0     0      1           ];

controlPts(:,1,2,2) = [ -L-R_i1   R_i1	 0    	1/sqrt(2)   ];
controlPts(:,2,2,2) = [ -L-R_i1   R_i1	 R_i1   	1/2         ];
controlPts(:,3,2,2) = [ -L-R_i1   0     R_i1 	1/sqrt(2)   ];
controlPts(:,4,2,2) = [ -L-R_i1  -R_i1   R_i1   	1/2         ];
controlPts(:,5,2,2) = [ -L-R_i1  -R_i1   0    	1/sqrt(2)   ];
controlPts(:,6,2,2) = [ -L-R_i1  -R_i1  -R_i1  	1/2         ];
controlPts(:,7,2,2) = [ -L-R_i1   0    -R_i1   	1/sqrt(2)   ];
controlPts(:,8,2,2) = [ -L-R_i1   R_i1  -R_i1   	1/2         ];
controlPts(:,9,2,2) = [ -L-R_i1   R_i1   0   	1/sqrt(2)   ];

controlPts(:,1,3,2) = [ -L     R_i1   0    	1           ];
controlPts(:,2,3,2) = [ -L     R_i1   R_i1   	1/sqrt(2)   ];
controlPts(:,3,3,2) = [	-L     0     R_i1  	1           ];
controlPts(:,4,3,2) = [	-L    -R_i1   R_i1  	1/sqrt(2)   ];
controlPts(:,5,3,2) = [	-L    -R_i1   0    	1           ];
controlPts(:,6,3,2) = [	-L    -R_i1  -R_i1   	1/sqrt(2)   ];
controlPts(:,7,3,2) = [	-L     0    -R_i1   	1           ];
controlPts(:,8,3,2) = [	-L     R_i1  -R_i1   	1/sqrt(2)   ];
controlPts(:,9,3,2) = [	-L     R_i1   0    	1           ];

controlPts(:,1,4,2) = [  -L/2  R_im   0    	1           ];
controlPts(:,2,4,2) = [  -L/2  R_im   R_im   	1/sqrt(2)   ];
controlPts(:,3,4,2) = [  -L/2  0     R_im  	1           ];
controlPts(:,4,4,2) = [  -L/2 -R_im   R_im  	1/sqrt(2)   ];
controlPts(:,5,4,2) = [  -L/2 -R_im   0    	1           ];
controlPts(:,6,4,2) = [  -L/2 -R_im  -R_im   	1/sqrt(2)   ];
controlPts(:,7,4,2) = [  -L/2  0    -R_im   	1           ];
controlPts(:,8,4,2) = [  -L/2  R_im  -R_im   	1/sqrt(2)   ];
controlPts(:,9,4,2) = [  -L/2  R_im   0    	1           ];

controlPts(:,1,5,2) = [  0     R_i2   0    	1           ];
controlPts(:,2,5,2) = [  0     R_i2   R_i2   	1/sqrt(2)   ];
controlPts(:,3,5,2) = [	 0     0     R_i2  	1           ];
controlPts(:,4,5,2) = [	 0    -R_i2   R_i2  	1/sqrt(2)   ];
controlPts(:,5,5,2) = [	 0    -R_i2   0    	1           ];
controlPts(:,6,5,2) = [	 0    -R_i2  -R_i2   	1/sqrt(2)   ];
controlPts(:,7,5,2) = [	 0     0    -R_i2   	1           ];
controlPts(:,8,5,2) = [	 0     R_i2  -R_i2   	1/sqrt(2)   ];
controlPts(:,9,5,2) = [	 0     R_i2   0    	1           ];

controlPts(:,1,6,2) = [  R_i2   R_i2 	 0   	1/sqrt(2)   ];
controlPts(:,2,6,2) = [  R_i2   R_i2 	 R_i2   	1/2         ];
controlPts(:,3,6,2) = [  R_i2   0     R_i2  	1/sqrt(2)   ];
controlPts(:,4,6,2) = [  R_i2  -R_i2 	 R_i2  	1/2         ];
controlPts(:,5,6,2) = [  R_i2  -R_i2 	 0    	1/sqrt(2)   ];
controlPts(:,6,6,2) = [  R_i2  -R_i2  -R_i2 	1/2         ];
controlPts(:,7,6,2) = [  R_i2   0    -R_i2  	1/sqrt(2)   ];
controlPts(:,8,6,2) = [  R_i2   R_i2  -R_i2  	1/2         ];
controlPts(:,9,6,2) = [  R_i2   R_i2   0  	1/sqrt(2)   ];

controlPts(:,1,7,2) = [  R_i2   0     0      1           ];
controlPts(:,2,7,2) = [  R_i2   0     0      1/sqrt(2)   ];
controlPts(:,3,7,2) = [  R_i2   0     0      1           ];
controlPts(:,4,7,2) = [  R_i2   0     0      1/sqrt(2)   ];
controlPts(:,5,7,2) = [  R_i2   0     0      1           ];
controlPts(:,6,7,2) = [  R_i2   0     0      1/sqrt(2)   ];
controlPts(:,7,7,2) = [  R_i2   0     0      1           ];
controlPts(:,8,7,2) = [  R_i2   0     0      1/sqrt(2)   ];
controlPts(:,9,7,2) = [  R_i2   0     0      1           ];


switch alignWithAxis
    case 'Xaxis' 
        % Nothing to be done
    case 'Yaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = temp;
    case 'Zaxis'
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = temp;
end

controlPts(1,:,:,:) = controlPts(1,:,:,:) + x_0(1);
controlPts(2,:,:,:) = controlPts(2,:,:,:) + x_0(2);
controlPts(3,:,:,:) = controlPts(3,:,:,:) + x_0(3);

nurbs = createNURBSobject(controlPts,{Xi, Eta, Zeta});
