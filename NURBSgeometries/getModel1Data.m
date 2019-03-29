function solid = getModel1Data(R_i, R_o, L, eta1, eta2)
if nargin < 4
    eta1 = 0.327581905944480;
    eta2 = 0.8;
end
t = R_o-R_i;

if R_o/2 <= t
    error('Input must satisfy R_i/2 <= t!')
end

%% Parameters for outer surface
y_1 = R_o/20*(8-sqrt(29));
z_1 = R_o/10*(2+sqrt(29));

z_2 = z_1-y_1/z_1*(3*R_o/4-y_1);
theta_1 = acos(4*y_1/(3*R_o));
w_1 = cos(theta_1/2);

%% Parameters for inner surface
z_5 = sqrt(R_i^2-t^2);
z_6 = z_5 - t/z_5*(R_i-t);
y_4 = (R_o/2+R_i)/2;

a = 1+(2*z_5/(R_o-2*t))^2;
b =  -(2*z_5/(R_o-2*t))^2*R_o;
c =   (2*z_5/(R_o-2*t))^2*R_o^2/4 - y_4^2;
y_3 = (-b-sqrt(b^2-4*a*c))/(2*a);
z_3 = sqrt(y_4^2-y_3^2);


z_4 = z_3-y_3/z_3*(y_4-y_3);
theta_2 = acos(y_3/y_4);
w_2 = cos(theta_2/2);
theta_3 = acos(t/R_i);
w_3 = cos(theta_3/2);

x_1 = sqrt(R_i^2-t^2);
x_2 = x_1 - t/x_1*(R_i-t);

Xi = [0 0 0 1 1 2 2 3 3 3]/3;
Eta = [0 0 0 eta1 eta1 eta2 eta2 1 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,7,7,2);

%% inner surface
controlPts(:,1,1,1) = [ t-L   R_o/2     0      1           ];
controlPts(:,2,1,1) = [ t-L   R_o/2     0      w_2         ];
controlPts(:,3,1,1) = [ t-L   R_o/2     0      1           ];
controlPts(:,4,1,1) = [ t-L   R_o/2     0      1           ];
controlPts(:,5,1,1) = [ t-L   R_o/2     0      1           ];
controlPts(:,6,1,1) = [ t-L   R_o/2     0      w_2         ];
controlPts(:,7,1,1) = [ t-L   R_o/2     0      1           ];

controlPts(:,1,2,1) = [ t-L     y_4      0          1           ];
controlPts(:,2,2,1) = [ t-L     y_4      z_4   	   w_2         ];
controlPts(:,3,2,1) = [	t-L     y_3      z_3  	   1           ];
controlPts(:,4,2,1) = [	t-L     y_3      0          1           ];
controlPts(:,5,2,1) = [	t-L     y_3     -z_3   	   1           ];
controlPts(:,6,2,1) = [ t-L     y_4     -z_4   	   w_2         ];
controlPts(:,7,2,1) = [ t-L     y_4      0          1           ];

controlPts(:,1,3,1) = [ t-L     R_i   0          1           ];
controlPts(:,2,3,1) = [ t-L     R_i   z_6        w_3         ];
controlPts(:,3,3,1) = [	t-L     t   	 z_5        1           ];
controlPts(:,4,3,1) = [	t-L     t     0          1           ];
controlPts(:,5,3,1) = [	t-L     t    -z_5        1           ];
controlPts(:,6,3,1) = [ t-L     R_i  -z_6        w_3         ];
controlPts(:,7,3,1) = [ t-L     R_i   0          1           ];

controlPts(:,1,4,1) = [  -L/2     R_i   0          1           ];
controlPts(:,2,4,1) = [  -L/2     R_i   z_6        w_3         ];
controlPts(:,3,4,1) = [	 -L/2     t   	 z_5        1           ];
controlPts(:,4,4,1) = [	 -L/2     t     0          1           ];
controlPts(:,5,4,1) = [	 -L/2     t    -z_5        1           ];
controlPts(:,6,4,1) = [  -L/2     R_i  -z_6        w_3         ];
controlPts(:,7,4,1) = [  -L/2     R_i   0          1           ];

controlPts(:,1,5,1) = [ 0     R_i   0          1           ];
controlPts(:,2,5,1) = [ 0     R_i   z_6        w_3         ];
controlPts(:,3,5,1) = [	0     t     z_5        1           ];
controlPts(:,4,5,1) = [	0     t     0          1           ];
controlPts(:,5,5,1) = [	0     t    -z_5        1           ];
controlPts(:,6,5,1) = [ 0     R_i  -z_6        w_3         ];
controlPts(:,7,5,1) = [ 0     R_i   0          1           ];

controlPts(:,1,6,1) = [ x_2     R_i   0          w_3           ];
controlPts(:,2,6,1) = [ x_2     R_i   z_6        w_3*w_3       ];
controlPts(:,3,6,1) = [	x_2     t     z_5        w_3           ];
controlPts(:,4,6,1) = [	x_2     t     0          w_3           ];
controlPts(:,5,6,1) = [	x_2     t    -z_5        w_3           ];
controlPts(:,6,6,1) = [ x_2     R_i  -z_6        w_3*w_3       ];
controlPts(:,7,6,1) = [ x_2     R_i   0          w_3           ];

controlPts(:,1,7,1) = [ x_1   t     0      1           ];
controlPts(:,2,7,1) = [ x_1   t     0      w_3         ];
controlPts(:,3,7,1) = [ x_1   t     0      1           ];
controlPts(:,4,7,1) = [ x_1   t     0      1           ];
controlPts(:,5,7,1) = [ x_1   t     0      1           ];
controlPts(:,6,7,1) = [ x_1   t     0      w_3         ];
controlPts(:,7,7,1) = [ x_1   t     0      1           ];


%% outer surface
controlPts(:,1,1,2) = [ -L   R_o/2     0      1           ];
controlPts(:,2,1,2) = [ -L   R_o/2     0      w_1         ];
controlPts(:,3,1,2) = [ -L   R_o/2     0      1           ];
controlPts(:,4,1,2) = [ -L   R_o/2     0      1           ];
controlPts(:,5,1,2) = [ -L   R_o/2     0      1           ];
controlPts(:,6,1,2) = [ -L   R_o/2     0      w_1         ];
controlPts(:,7,1,2) = [ -L   R_o/2     0      1           ];

controlPts(:,1,2,2) = [ -L     3*R_o/4   0        1           ];
controlPts(:,2,2,2) = [ -L     3*R_o/4   z_2   	w_1         ];
controlPts(:,3,2,2) = [	-L     y_1       z_1  	1           ];
controlPts(:,4,2,2) = [	-L     y_1       0        1           ];
controlPts(:,5,2,2) = [	-L     y_1      -z_1   	1           ];
controlPts(:,6,2,2) = [ -L     3*R_o/4  -z_2   	w_1         ];
controlPts(:,7,2,2) = [ -L     3*R_o/4   0        1           ];

controlPts(:,1,3,2) = [ -L     R_o   0            1           ];
controlPts(:,2,3,2) = [ -L     R_o   R_o          1/sqrt(2)   ];
controlPts(:,3,3,2) = [	-L     0     R_o          1           ];
controlPts(:,4,3,2) = [	-L     0     0            1           ];
controlPts(:,5,3,2) = [	-L     0    -R_o         	1           ];
controlPts(:,6,3,2) = [ -L     R_o  -R_o          1/sqrt(2)   ];
controlPts(:,7,3,2) = [ -L     R_o   0            1           ];

controlPts(:,1,4,2) = [ -L/2     R_o   0       1              ];
controlPts(:,2,4,2) = [ -L/2     R_o   R_o   	1/sqrt(2)       ];
controlPts(:,3,4,2) = [	-L/2     0     R_o  	1               ];
controlPts(:,4,4,2) = [	-L/2     0     0       1              ];
controlPts(:,5,4,2) = [	-L/2     0    -R_o   	1               ];
controlPts(:,6,4,2) = [ -L/2     R_o  -R_o   	1/sqrt(2)       ];
controlPts(:,7,4,2) = [ -L/2     R_o   0       1              ];

controlPts(:,1,5,2) = [ 0     R_o   0         1           ];
controlPts(:,2,5,2) = [ 0     R_o   R_o   	1/sqrt(2)   ];
controlPts(:,3,5,2) = [	0     0     R_o       1           ];
controlPts(:,4,5,2) = [	0     0     0         1           ];
controlPts(:,5,5,2) = [	0     0    -R_o   	1           ];
controlPts(:,6,5,2) = [ 0     R_o  -R_o       1/sqrt(2)   ];
controlPts(:,7,5,2) = [ 0     R_o   0         1           ];

controlPts(:,1,6,2) = [ R_o     R_o   0         1/sqrt(2)	];
controlPts(:,2,6,2) = [ R_o     R_o   R_o   	1/2      	];
controlPts(:,3,6,2) = [	R_o     0     R_o  	    1/sqrt(2)	];
controlPts(:,4,6,2) = [	R_o     0     0         1	];
controlPts(:,5,6,2) = [	R_o     0    -R_o   	1/sqrt(2)   ];
controlPts(:,6,6,2) = [ R_o     R_o  -R_o   	1/2      	];
controlPts(:,7,6,2) = [ R_o     R_o   0         1/sqrt(2) 	];

controlPts(:,1,7,2) = [ R_o   0     0      1           ];
controlPts(:,2,7,2) = [ R_o   0     0      1/sqrt(2)   ];
controlPts(:,3,7,2) = [ R_o   0     0      1           ];
controlPts(:,4,7,2) = [ R_o   0     0      1           ];
controlPts(:,5,7,2) = [ R_o   0     0      1           ];
controlPts(:,6,7,2) = [ R_o   0     0      1/sqrt(2)   ];
controlPts(:,7,7,2) = [ R_o   0     0      1           ];

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
