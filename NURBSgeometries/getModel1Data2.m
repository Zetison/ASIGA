function solid = getModel1Data2(R_o,L, x_0, alignWithAxis, shift, xi1, eta1, eta2)

if nargin < 5
    alignWithAxis = 'Xaxis';
end

%% Parameters for outer surface
y_1 = R_o/20*(8-sqrt(29));
z_1 = R_o/10*(2+sqrt(29));

z_2 = z_1-y_1/z_1*(3*R_o/4-y_1);
theta_1 = acos(4*y_1/(3*R_o));
w_1 = cos(theta_1/2);


Xi = [0 0 0 1 1 1+xi1 1+xi1 2-xi1 2-xi1 2 2 3 3 3]/3;
Eta = [0 0 0 eta1 eta1 eta2 eta2 1 1 1];

controlPts = zeros(4,11,7);

%% outer surface
controlPts(:,1,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,2,1) = [ -L-shift   R_o/2     0      w_1         ];
controlPts(:,3,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,4,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,5,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,6,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,7,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,8,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,9,1) = [ -L-shift   R_o/2     0      1           ];
controlPts(:,10,1) = [ -L-shift   R_o/2     0      w_1         ];
controlPts(:,11,1) = [ -L-shift   R_o/2     0      1           ];

controlPts(:,1,2) = [ -L-shift     3*R_o/4   0        1           ];
controlPts(:,2,2) = [ -L-shift     3*R_o/4   z_2   	w_1         ];
controlPts(:,3,2) = [	-L-shift     y_1       z_1  	1           ];
controlPts(:,4,2) = [	-L-shift     y_1-shift/4       z_1  	1           ];
controlPts(:,5,2) = [	-L-shift     y_1-shift/2       z_1  	1           ];
controlPts(:,6,2) = [	-L-shift     y_1-shift/2       0        1           ];
controlPts(:,7,2) = [	-L-shift     y_1-shift/2      -z_1   	1           ];
controlPts(:,8,2) = [	-L-shift     y_1-shift/4      -z_1   	1           ];
controlPts(:,9,2) = [	-L-shift     y_1      -z_1   	1           ];
controlPts(:,10,2) = [ -L-shift     3*R_o/4  -z_2   	w_1         ];
controlPts(:,11,2) = [ -L-shift     3*R_o/4   0        1           ];

controlPts(:,1,3) = [ -L-shift     R_o   0            1           ];
controlPts(:,2,3) = [ -L-shift     R_o   R_o          1/sqrt(2)   ];
controlPts(:,3,3) = [	-L-shift     0     R_o          1           ];
controlPts(:,4,3) = [	-L-shift     -shift/2     R_o          1           ];
controlPts(:,5,3) = [	-L-shift     -shift     R_o          1           ];
controlPts(:,6,3) = [	-L-shift     -shift     0            1           ];
controlPts(:,7,3) = [	-L-shift     -shift    -R_o         	1           ];
controlPts(:,8,3) = [	-L-shift     -shift/2    -R_o         	1           ];
controlPts(:,9,3) = [	-L-shift     0    -R_o         	1           ];
controlPts(:,10,3) = [ -L-shift     R_o  -R_o          1/sqrt(2)   ];
controlPts(:,11,3) = [ -L-shift     R_o   0            1           ];

controlPts(:,1,4) = [ -L/2     R_o   0       1              ];
controlPts(:,2,4) = [ -L/2     R_o   R_o   	1/sqrt(2)       ];
controlPts(:,3,4) = [	-L/2     0     R_o  	1               ];
controlPts(:,4,4) = [	-L/2     -shift/2     R_o  	1               ];
controlPts(:,5,4) = [	-L/2     -shift     R_o  	1               ];
controlPts(:,6,4) = [	-L/2     -shift     0       1              ];
controlPts(:,7,4) = [	-L/2     -shift    -R_o   	1               ];
controlPts(:,8,4) = [	-L/2     -shift/2    -R_o   	1               ];
controlPts(:,9,4) = [	-L/2     0    -R_o   	1               ];
controlPts(:,10,4) = [ -L/2     R_o  -R_o   	1/sqrt(2)       ];
controlPts(:,11,4) = [ -L/2     R_o   0       1              ];

controlPts(:,1,5) = [ 0     R_o   0         1           ];
controlPts(:,2,5) = [ 0     R_o   R_o   	1/sqrt(2)   ];
controlPts(:,3,5) = [	0     0     R_o       1           ];
controlPts(:,4,5) = [	0     -shift/2     R_o       1           ];
controlPts(:,5,5) = [	0     -shift     R_o       1           ];
controlPts(:,6,5) = [	0     -shift     0         1           ];
controlPts(:,7,5) = [	0     -shift    -R_o   	1           ];
controlPts(:,8,5) = [	0     -shift/2    -R_o   	1           ];
controlPts(:,9,5) = [	0     0    -R_o   	1           ];
controlPts(:,10,5) = [ 0     R_o  -R_o       1/sqrt(2)   ];
controlPts(:,11,5) = [ 0     R_o   0         1           ];

controlPts(:,1,6) = [ R_o     R_o   0         1/sqrt(2)	];
controlPts(:,2,6) = [ R_o     R_o   R_o   	1/2      	];
controlPts(:,3,6) = [	R_o     0     R_o  	    1/sqrt(2)	];
controlPts(:,4,6) = [	R_o     -shift/2     R_o       1/sqrt(2)           ];
controlPts(:,5,6) = [	R_o     -shift     R_o       1/sqrt(2)           ];
controlPts(:,6,6) = [	R_o     -shift     0         1/sqrt(2)           ];
controlPts(:,7,6) = [	R_o     -shift    -R_o   	1/sqrt(2)           ];
controlPts(:,8,6) = [	R_o     -shift/2    -R_o   	1/sqrt(2)           ];
controlPts(:,9,6) = [	R_o     0    -R_o   	1/sqrt(2)   ];
controlPts(:,10,6) = [ R_o     R_o  -R_o   	1/2      	];
controlPts(:,11,6) = [ R_o     R_o   0         1/sqrt(2) 	];

controlPts(:,1,7) = [ R_o   0     0      1           ];
controlPts(:,2,7) = [ R_o   0     0      1/sqrt(2)   ];
controlPts(:,3,7) = [ R_o   0     0      1           ];
controlPts(:,4,7) = [ R_o   -shift/2     0      1           ];
controlPts(:,5,7) = [ R_o   -shift     0      1           ];
controlPts(:,6,7) = [ R_o   -shift     0      1           ];
controlPts(:,7,7) = [ R_o   -shift     0      1           ];
controlPts(:,8,7) = [ R_o   -shift/2     0      1           ];
controlPts(:,9,7) = [ R_o   0     0      1           ];
controlPts(:,10,7) = [ R_o   0     0      1/sqrt(2)   ];
controlPts(:,11,7) = [ R_o   0     0      1           ];

% controlPts(2,:,:,:) = controlPts(2,:,:,:) - R_o/2;

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

solid = createNURBSobject(controlPts,{Xi, Eta});
