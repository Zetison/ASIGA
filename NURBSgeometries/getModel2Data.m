function solid = getModel2Data(R_o,R_i,L)

t = R_o-R_i;

x1 = L/4;

w = cos(pi/8);

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Zeta = [0 0 1 1];

controlPts = zeros(4,9,9,2);

% inner surface
% controlPts(:,1,1,1) = [ -L/2+t   0     0      1           ];
% controlPts(:,2,1,1) = [ -L/2+t   0     0      1/sqrt(2)   ];
% controlPts(:,3,1,1) = [ -L/2+t   0     0      1           ];
% controlPts(:,4,1,1) = [ -L/2+t   0     0      1/sqrt(2)   ];
% controlPts(:,5,1,1) = [ -L/2+t   0     0      1           ];
% controlPts(:,6,1,1) = [ -L/2+t   0     0      1/sqrt(2)   ];
% controlPts(:,7,1,1) = [ -L/2+t   0     0      1           ];
% controlPts(:,8,1,1) = [ -L/2+t   0     0      1/sqrt(2)   ];
% controlPts(:,9,1,1) = [ -L/2+t   0     0      1           ];
% 
% controlPts(:,1,2,1) = [ -L/2+t   R_i/2	 0    	1   ];
% controlPts(:,2,2,1) = [ -L/2+t   R_i/2	 R_i/2   	1/sqrt(2)         ];
% controlPts(:,3,2,1) = [ -L/2+t   0     R_i/2 	1   ];
% controlPts(:,4,2,1) = [ -L/2+t  -R_i/2   R_i/2   	1/sqrt(2)        ];
% controlPts(:,5,2,1) = [ -L/2+t  -R_i/2   0    	1   ];
% controlPts(:,6,2,1) = [ -L/2+t  -R_i/2  -R_i/2  	1/sqrt(2)        ];
% controlPts(:,7,2,1) = [ -L/2+t   0    -R_i/2   	1   ];
% controlPts(:,8,2,1) = [ -L/2+t   R_i/2  -R_i/2   	1/sqrt(2)         ];
% controlPts(:,9,2,1) = [ -L/2+t   R_i/2   0   	1   ];
% 
% controlPts(:,1,3,1) = [ -L/2+t     R_i   0    	1           ];
% controlPts(:,2,3,1) = [ -L/2+t     R_i   R_i   	1/sqrt(2)   ];
% controlPts(:,3,3,1) = [	-L/2+t     0     R_i  	1           ];
% controlPts(:,4,3,1) = [	-L/2+t    -R_i   R_i  	1/sqrt(2)   ];
% controlPts(:,5,3,1) = [	-L/2+t    -R_i   0    	1           ];
% controlPts(:,6,3,1) = [	-L/2+t    -R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,7,3,1) = [	-L/2+t     0    -R_i   	1           ];
% controlPts(:,8,3,1) = [	-L/2+t     R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,9,3,1) = [	-L/2+t     R_i   0    	1           ];
% 
% controlPts(:,1,4,1) = [  0     R_i   0    	1           ];
% controlPts(:,2,4,1) = [  0     R_i   R_i   	1/sqrt(2)   ];
% controlPts(:,3,4,1) = [	 0     0     R_i  	1           ];
% controlPts(:,4,4,1) = [	 0    -R_i   R_i  	1/sqrt(2)   ];
% controlPts(:,5,4,1) = [	 0    -R_i   0    	1           ];
% controlPts(:,6,4,1) = [	 0    -R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,7,4,1) = [	 0     0    -R_i   	1           ];
% controlPts(:,8,4,1) = [	 0     R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,9,4,1) = [	 0     R_i   0    	1           ];
% 
% controlPts(:,1,5,1) = [  L/2     R_i   0    	1           ];
% controlPts(:,2,5,1) = [  L/2     R_i   R_i   	1/sqrt(2)   ];
% controlPts(:,3,5,1) = [	 L/2     0     R_i  	1           ];
% controlPts(:,4,5,1) = [	 L/2    -R_i   R_i  	1/sqrt(2)   ];
% controlPts(:,5,5,1) = [	 L/2    -R_i   0    	1           ];
% controlPts(:,6,5,1) = [	 L/2    -R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,7,5,1) = [	 L/2     0    -R_i   	1           ];
% controlPts(:,8,5,1) = [	 L/2     R_i  -R_i   	1/sqrt(2)   ];
% controlPts(:,9,5,1) = [	 L/2     R_i   0    	1           ];
% 
% controlPts(:,1,6,1) = [  L/2+R_i   R_i 	 0   	1/sqrt(2)   ];
% controlPts(:,2,6,1) = [  L/2+R_i   R_i 	 R_i   	1/2         ];
% controlPts(:,3,6,1) = [  L/2+R_i   0     R_i  	1/sqrt(2)   ];
% controlPts(:,4,6,1) = [  L/2+R_i  -R_i 	 R_i  	1/2         ];
% controlPts(:,5,6,1) = [  L/2+R_i  -R_i 	 0    	1/sqrt(2)   ];
% controlPts(:,6,6,1) = [  L/2+R_i  -R_i  -R_i 	1/2         ];
% controlPts(:,7,6,1) = [  L/2+R_i   0    -R_i  	1/sqrt(2)   ];
% controlPts(:,8,6,1) = [  L/2+R_i   R_i  -R_i  	1/2         ];
% controlPts(:,9,6,1) = [  L/2+R_i   R_i   0  	1/sqrt(2)   ];
% 
% controlPts(:,1,7,1) = [  L/2+R_i   0     0      1           ];
% controlPts(:,2,7,1) = [  L/2+R_i   0     0      1/sqrt(2)   ];
% controlPts(:,3,7,1) = [  L/2+R_i   0     0      1           ];
% controlPts(:,4,7,1) = [  L/2+R_i   0     0      1/sqrt(2)   ];
% controlPts(:,5,7,1) = [  L/2+R_i   0     0      1           ];
% controlPts(:,6,7,1) = [  L/2+R_i   0     0      1/sqrt(2)   ];
% controlPts(:,7,7,1) = [  L/2+R_i   0     0      1           ];
% controlPts(:,8,7,1) = [  L/2+R_i   0     0      1/sqrt(2)   ];
% controlPts(:,9,7,1) = [  L/2+R_i   0     0      1           ];

% outer surface
controlPts(:,1,1,2) = [ -L/2   0     0      1           ];
controlPts(:,2,1,2) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts(:,3,1,2) = [ -L/2   0     0      1           ];
controlPts(:,4,1,2) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts(:,5,1,2) = [ -L/2   0     0      1           ];
controlPts(:,6,1,2) = [ -L/2   0     0      w           ];
controlPts(:,7,1,2) = [ -L/2   0     0      1           ];
controlPts(:,8,1,2) = [ -L/2   0     0      w           ];
controlPts(:,9,1,2) = [ -L/2   0     0      1           ];

controlPts(:,1,2,2) = [ -L/2   R_o/2	 0    	1   ];
controlPts(:,2,2,2) = [ -L/2   R_o/2	 R_o/2  1/sqrt(2)         ];
controlPts(:,3,2,2) = [ -L/2   0        R_o/2 	1   ];
controlPts(:,4,2,2) = [ -L/2  -R_o/2   R_o/2   	1/sqrt(2)        ];
controlPts(:,5,2,2) = [ -L/2  -R_o/2   0    	1   ];
controlPts(:,6,2,2) = [ -L/2  -R_o/2  -R_o/2  	1/sqrt(2)         ];
controlPts(:,7,2,2) = [ -L/2   0    -R_o/2   	1   ];
controlPts(:,8,2,2) = [ -L/2   R_o/2  -R_o/2   	1/sqrt(2)         ];
controlPts(:,9,2,2) = [ -L/2   R_o/2   0        1   ];

controlPts(:,1,3,2) = [ -L/2     R_o   0    	1           ];
controlPts(:,2,3,2) = [ -L/2     R_o   R_o   	1/sqrt(2)   ];
controlPts(:,3,3,2) = [	-L/2     0     R_o  	1           ];
controlPts(:,4,3,2) = [	-L/2    -R_o   R_o  	1/sqrt(2)   ];
controlPts(:,5,3,2) = [	-L/2    -R_o   0    	1           ];
controlPts(:,6,3,2) = [	-L/2    -R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,7,3,2) = [	-L/2     0    -R_o   	1           ];
controlPts(:,8,3,2) = [	-L/2     R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,9,3,2) = [	-L/2     R_o   0    	1           ];

controlPts(:,1,4,2) = [ 0     R_o   0    	1           ];
controlPts(:,2,4,2) = [ 0     R_o   R_o   	1/sqrt(2)   ];
controlPts(:,3,4,2) = [	0     0     R_o  	1           ];
controlPts(:,4,4,2) = [	0    -R_o   R_o  	1/sqrt(2)   ];
controlPts(:,5,4,2) = [	0    -R_o   0    	1           ];
controlPts(:,6,4,2) = [	0    -R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,7,4,2) = [	0     0    -R_o   	1           ];
controlPts(:,8,4,2) = [	0     R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,9,4,2) = [	0     R_o   0    	1           ];

controlPts(:,1,5,2) = [  L/2     R_o   0    	1           ];
controlPts(:,2,5,2) = [  L/2     R_o   R_o   	1/sqrt(2)   ];
controlPts(:,3,5,2) = [	 L/2     0     R_o  	1           ];
controlPts(:,4,5,2) = [	 x1    -R_o   R_o  	1/sqrt(2)   ];
controlPts(:,5,5,2) = [	 L/2   -R_o   0    	1           ];
controlPts(:,6,5,2) = [	 L/2    -R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,7,5,2) = [	 L/2     0     -R_o   	1           ];
controlPts(:,8,5,2) = [	 L/2     R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,9,5,2) = [	 L/2     R_o   0    	1           ];

controlPts(:,1,6,2) = [  L/2+R_o         R_o   0    	1/sqrt(2)          ];
controlPts(:,2,6,2) = [  L/2+R_o         R_o   R_o   	1/2   ];
controlPts(:,3,6,2) = [	 L/2             0     R_o  	1/sqrt(2)           ];
controlPts(:,4,6,2) = [	 x1/2    -R_o   R_o  	1/sqrt(2)   ];
controlPts(:,5,6,2) = [	 L/2     -R_o   0    	1           ];
controlPts(:,6,6,2) = [	 L/2     -R_o  -R_o   	1/sqrt(2)   ];
controlPts(:,7,6,2) = [	 L/2             0    -R_o   	1/sqrt(2)           ];
controlPts(:,8,6,2) = [	 L/2+R_o        R_o  -R_o   	1/2   ];
controlPts(:,9,6,2) = [	 L/2+R_o        R_o   0    	1/sqrt(2)           ];

controlPts(:,1,7,2) = [  L/2+R_o                 0               	0       1           ];
controlPts(:,2,7,2) = [  L/2+R_o                 0              	R_o     1/sqrt(2)   ];
controlPts(:,3,7,2) = [  L/2                     0                  R_o 	1           ];
controlPts(:,4,7,2) = [  L/2                    -R_o                R_o  	1/sqrt(2)   ];
controlPts(:,5,7,2) = [  L/2                    -R_o                0       1           ];
controlPts(:,6,7,2) = [  L/2+R_o*(sqrt(2)-1)	-R_o                0       w           ];
controlPts(:,7,7,2) = [  L/2+R_o/sqrt(2)        -R_o/sqrt(2)        0       1           ];
controlPts(:,8,7,2) = [  L/2+R_o                 R_o*(1-sqrt(2))	0       w           ];
controlPts(:,9,7,2) = [  L/2+R_o                 0                 	0       1           ];

controlPts(:,1,8,2) = [  L/2+R_o/2                  0                   0       1           ];
controlPts(:,2,8,2) = [  L/2+R_o/2                  0                   R_o/2	1/sqrt(2)   ];
controlPts(:,3,8,2) = [  L/2                        0                   R_o/2 	1           ];
controlPts(:,4,8,2) = [  L/2                        -R_o/2              R_o/2  	1/sqrt(2)   ];
controlPts(:,5,8,2) = [  L/2                        -R_o/2              0       1           ];
controlPts(:,6,8,2) = [  L/2+R_o*(1/sqrt(2)-1/2)	-R_o/2              0       w           ];
controlPts(:,7,8,2) = [  L/2+R_o/(2*sqrt(2))        -R_o/(2*sqrt(2)) 	0       1           ];
controlPts(:,8,8,2) = [  L/2+R_o/2                  R_o*(1/2-1/sqrt(2))	0       w           ];
controlPts(:,9,8,2) = [  L/2+R_o/2                  0                   0       1           ];

controlPts(:,1,9,2) = [  L/2+R_o   0     0      1           ];
controlPts(:,2,9,2) = [  L/2+R_o   0     0      1/sqrt(2)   ];
controlPts(:,3,9,2) = [  L/2+R_o   0     0      1           ];
controlPts(:,4,9,2) = [  L/2+R_o   0     0      w   ];
controlPts(:,5,9,2) = [  L/2+R_o   0     0      1           ];
controlPts(:,6,9,2) = [  L/2+R_o   0     0      w           ];
controlPts(:,7,9,2) = [  L/2+R_o   0     0      1           ];
controlPts(:,8,9,2) = [  L/2+R_o   0     0      1/sqrt(2)           ];
controlPts(:,9,9,2) = [  L/2+R_o   0     0      1           ];
controlPts(:,1,9,2) = [  L/2+R_o   0     0      1/sqrt(2)           ];
controlPts(:,2,9,2) = [  L/2+R_o   0     0      1   ];



solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
