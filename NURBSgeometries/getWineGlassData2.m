function nurbs = getWineGlassData2()

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0 0 0 0.5 1 1 1.1 1.9 2 2 2.2 2.5 3.2 4 4 4]/4;
Zeta = [0 0 1 1];

t1 = 0.001;
t2 = 0.0035; % Thickness of base

k = 0.001;

h1 = 0.005;
h2 = 0.115;
h6 = 0.125;
h7 = 0.15; 
h8 = 0.2;
h9 = 0.22;

r = 0.017/2; % Shaft radius at top and bottom

x_1 = 0.6*r;
y_1 = h1+0.6*r;

x_2 = x_1;
y_2 = h2-0.6*r;


R_0i = 0.073/2; % messured
R_1i = 0.02;
R_2i = 0.002;
R_5i = 0.005;   % bottom of cup
R_8i = 0.074/2; % messured
R_6i = R_8i - 0.0;
R_7i = R_8i + 0.017;

R_0o = R_0i;
R_1o = 0.6*R_0o;
R_4o = 0.2*x_1; %bending radius of shaft
R_6o = R_6i+t1;
R_7o = R_7i+t1;
R_8o = R_8i+t1;

controlPts = zeros(4,9,13,2);

controlPts(:,1,1,1) = [  R_0i	0       0	 1           ];
controlPts(:,2,1,1) = [  R_0i	R_0i    0	 1/sqrt(2)   ];
controlPts(:,3,1,1) = [  0      R_0i    0    1           ];
controlPts(:,4,1,1) = [ -R_0i	R_0i    0    1/sqrt(2)   ];
controlPts(:,5,1,1) = [ -R_0i   0       0    1           ];
controlPts(:,6,1,1) = [ -R_0i  -R_0i    0    1/sqrt(2)   ];
controlPts(:,7,1,1) = [  0     -R_0i    0    1           ];
controlPts(:,8,1,1) = [  R_0i  -R_0i    0    1/sqrt(2)   ];
controlPts(:,9,1,1) = [  R_0i	0       0    1           ];

controlPts(:,1,2,1) = [  R_1i	0       0	 1           ];
controlPts(:,2,2,1) = [  R_1i	R_1i    0	 1/sqrt(2)   ];
controlPts(:,3,2,1) = [  0      R_1i    0    1           ];
controlPts(:,4,2,1) = [ -R_1i	R_1i    0    1/sqrt(2)   ];
controlPts(:,5,2,1) = [ -R_1i   0       0    1           ];
controlPts(:,6,2,1) = [ -R_1i  -R_1i    0    1/sqrt(2)   ];
controlPts(:,7,2,1) = [  0     -R_1i    0    1           ];
controlPts(:,8,2,1) = [  R_1i  -R_1i    0    1/sqrt(2)   ];
controlPts(:,9,2,1) = [  R_1i	0       0    1           ];

controlPts(:,1,3,1) = [  R_2i	0       h1	  1           ];
controlPts(:,2,3,1) = [  R_2i	R_2i    h1	  1/sqrt(2)   ];
controlPts(:,3,3,1) = [  0      R_2i    h1    1           ];
controlPts(:,4,3,1) = [ -R_2i	R_2i    h1    1/sqrt(2)   ];
controlPts(:,5,3,1) = [ -R_2i   0       h1    1           ];
controlPts(:,6,3,1) = [ -R_2i  -R_2i    h1    1/sqrt(2)   ];
controlPts(:,7,3,1) = [  0     -R_2i    h1    1           ];
controlPts(:,8,3,1) = [  R_2i  -R_2i    h1    1/sqrt(2)   ];
controlPts(:,9,3,1) = [  R_2i	0       h1    1           ];

controlPts(:,1,4,1) = [  0	  0     h1    1           ];
controlPts(:,2,4,1) = [  0	  0     h1    1/sqrt(2)   ];
controlPts(:,3,4,1) = [  0    0     h1    1           ];
controlPts(:,4,4,1) = [  0    0     h1    1/sqrt(2)   ];
controlPts(:,5,4,1) = [  0    0     h1    1           ];
controlPts(:,6,4,1) = [  0    0     h1    1/sqrt(2)   ];
controlPts(:,7,4,1) = [  0    0     h1    1           ];
controlPts(:,8,4,1) = [  0    0     h1    1/sqrt(2)   ];
controlPts(:,9,4,1) = [  0    0     h1    1           ];

controlPts(:,1,5,1) = [  0	  0     y_1+k    1           ];
controlPts(:,2,5,1) = [  0	  0     y_1+k    1/sqrt(2)   ];
controlPts(:,3,5,1) = [  0    0     y_1+k    1           ];
controlPts(:,4,5,1) = [  0    0     y_1+k    1/sqrt(2)   ];
controlPts(:,5,5,1) = [  0    0     y_1+k    1           ];
controlPts(:,6,5,1) = [  0    0     y_1+k    1/sqrt(2)   ];
controlPts(:,7,5,1) = [  0    0     y_1+k    1           ];
controlPts(:,8,5,1) = [  0    0     y_1+k    1/sqrt(2)   ];
controlPts(:,9,5,1) = [  0    0     y_1+k    1           ];

controlPts(:,1,6,1) = [  0	  0     2*(h1+h2)/4    1           ];
controlPts(:,2,6,1) = [  0	  0     2*(h1+h2)/4    1/sqrt(2)   ];
controlPts(:,3,6,1) = [  0    0     2*(h1+h2)/4    1           ];
controlPts(:,4,6,1) = [  0    0     2*(h1+h2)/4    1/sqrt(2)   ];
controlPts(:,5,6,1) = [  0    0     2*(h1+h2)/4    1           ];
controlPts(:,6,6,1) = [  0    0     2*(h1+h2)/4    1/sqrt(2)   ];
controlPts(:,7,6,1) = [  0    0     2*(h1+h2)/4    1           ];
controlPts(:,8,6,1) = [  0    0     2*(h1+h2)/4    1/sqrt(2)   ];
controlPts(:,9,6,1) = [  0    0     2*(h1+h2)/4    1           ];

controlPts(:,1,7,1) = [  0	  0     y_2-k    1           ];
controlPts(:,2,7,1) = [  0	  0     y_2-k    1/sqrt(2)   ];
controlPts(:,3,7,1) = [  0    0     y_2-k    1           ];
controlPts(:,4,7,1) = [  0    0     y_2-k    1/sqrt(2)   ];
controlPts(:,5,7,1) = [  0    0     y_2-k    1           ];
controlPts(:,6,7,1) = [  0    0     y_2-k    1/sqrt(2)   ];
controlPts(:,7,7,1) = [  0    0     y_2-k    1           ];
controlPts(:,8,7,1) = [  0    0     y_2-k    1/sqrt(2)   ];
controlPts(:,9,7,1) = [  0    0     y_2-k    1           ];

controlPts(:,1,8,1) = [  0	  0     h2    1           ];
controlPts(:,2,8,1) = [  0	  0     h2    1/sqrt(2)   ];
controlPts(:,3,8,1) = [  0    0     h2    1           ];
controlPts(:,4,8,1) = [  0    0     h2    1/sqrt(2)   ];
controlPts(:,5,8,1) = [  0    0     h2    1           ];
controlPts(:,6,8,1) = [  0    0     h2    1/sqrt(2)   ];
controlPts(:,7,8,1) = [  0    0     h2    1           ];
controlPts(:,8,8,1) = [  0    0     h2    1/sqrt(2)   ];
controlPts(:,9,8,1) = [  0    0     h2    1           ];

controlPts(:,1,9,1) = [  R_5i	0       h2    1           ];
controlPts(:,2,9,1) = [  R_5i	R_5i    h2    1/sqrt(2)   ];
controlPts(:,3,9,1) = [  0  	R_5i    h2    1           ];
controlPts(:,4,9,1) = [ -R_5i    R_5i   h2    1/sqrt(2)   ];
controlPts(:,5,9,1) = [ -R_5i    0      h2    1           ];
controlPts(:,6,9,1) = [ -R_5i   -R_5i   h2    1/sqrt(2)   ];
controlPts(:,7,9,1) = [  0     -R_5i    h2    1           ];
controlPts(:,8,9,1) = [  R_5i   -R_5i   h2    1/sqrt(2)   ];
controlPts(:,9,9,1) = [  R_5i  	0       h2    1           ];

controlPts(:,1,10,1) = [  R_6i	0       h6+t1/2    1           ];
controlPts(:,2,10,1) = [  R_6i	R_6i    h6+t1/2    1/sqrt(2)   ];
controlPts(:,3,10,1) = [  0  	R_6i    h6+t1/2    1           ];
controlPts(:,4,10,1) = [ -R_6i    R_6i  h6+t1/2    1/sqrt(2)   ];
controlPts(:,5,10,1) = [ -R_6i    0     h6+t1/2    1           ];
controlPts(:,6,10,1) = [ -R_6i   -R_6i  h6+t1/2    1/sqrt(2)   ];
controlPts(:,7,10,1) = [  0     -R_6i   h6+t1/2    1           ];
controlPts(:,8,10,1) = [  R_6i   -R_6i  h6+t1/2    1/sqrt(2)   ];
controlPts(:,9,10,1) = [  R_6i  	0   h6+t1/2    1           ];

controlPts(:,1,11,1) = [  R_7i	0       h7    1           ];
controlPts(:,2,11,1) = [  R_7i	R_7i    h7    1/sqrt(2)   ];
controlPts(:,3,11,1) = [  0  	R_7i    h7    1           ];
controlPts(:,4,11,1) = [ -R_7i    R_7i   h7    1/sqrt(2)   ];
controlPts(:,5,11,1) = [ -R_7i    0      h7    1           ];
controlPts(:,6,11,1) = [ -R_7i   -R_7i   h7    1/sqrt(2)   ];
controlPts(:,7,11,1) = [  0     -R_7i    h7    1           ];
controlPts(:,8,11,1) = [  R_7i   -R_7i   h7    1/sqrt(2)   ];
controlPts(:,9,11,1) = [  R_7i  	0       h7    1           ];

controlPts(:,1,12,1) = [  R_8i	0        h8    1           ];
controlPts(:,2,12,1) = [  R_8i	R_8i     h8    1/sqrt(2)   ];
controlPts(:,3,12,1) = [  0  	R_8i     h8   1           ];
controlPts(:,4,12,1) = [ -R_8i    R_8i   h8   1/sqrt(2)   ];
controlPts(:,5,12,1) = [ -R_8i    0      h8    1           ];
controlPts(:,6,12,1) = [ -R_8i   -R_8i   h8    1/sqrt(2)   ];
controlPts(:,7,12,1) = [  0     -R_8i    h8    1           ];
controlPts(:,8,12,1) = [  R_8i   -R_8i   h8    1/sqrt(2)   ];
controlPts(:,9,12,1) = [  R_8i  	0    h8    1           ];

controlPts(:,1,13,1) = [  R_8i	0        h9    1           ];
controlPts(:,2,13,1) = [  R_8i	R_8i     h9    1/sqrt(2)   ];
controlPts(:,3,13,1) = [  0  	R_8i     h9   1           ];
controlPts(:,4,13,1) = [ -R_8i    R_8i   h9   1/sqrt(2)   ];
controlPts(:,5,13,1) = [ -R_8i    0      h9    1           ];
controlPts(:,6,13,1) = [ -R_8i   -R_8i   h9    1/sqrt(2)   ];
controlPts(:,7,13,1) = [  0     -R_8i    h9    1           ];
controlPts(:,8,13,1) = [  R_8i   -R_8i   h9    1/sqrt(2)   ];
controlPts(:,9,13,1) = [  R_8i  	0    h9    1           ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlPts(:,1,1,2) = [  R_0o	0        t2	  1           ];
controlPts(:,2,1,2) = [  R_0o	R_0o     t2	  1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0      R_0o     t2    1           ];
controlPts(:,4,1,2) = [ -R_0o   R_0o     t2    1/sqrt(2)   ];
controlPts(:,5,1,2) = [ -R_0o   0        t2    1           ];
controlPts(:,6,1,2) = [ -R_0o  -R_0o     t2    1/sqrt(2)   ];
controlPts(:,7,1,2) = [  0     -R_0o     t2    1           ];
controlPts(:,8,1,2) = [  R_0o  -R_0o     t2    1/sqrt(2)   ];
controlPts(:,9,1,2) = [  R_0o	0        t2    1           ];

controlPts(:,1,2,2) = [  R_1o	0        t2	  1           ];
controlPts(:,2,2,2) = [  R_1o	R_1o     t2	  1/sqrt(2)   ];
controlPts(:,3,2,2) = [  0      R_1o     t2    1           ];
controlPts(:,4,2,2) = [ -R_1o   R_1o     t2    1/sqrt(2)   ];
controlPts(:,5,2,2) = [ -R_1o   0        t2    1           ];
controlPts(:,6,2,2) = [ -R_1o  -R_1o     t2    1/sqrt(2)   ];
controlPts(:,7,2,2) = [  0     -R_1o     t2    1           ];
controlPts(:,8,2,2) = [  R_1o  -R_1o     t2    1/sqrt(2)   ];
controlPts(:,9,2,2) = [  R_1o	0        t2    1           ];

controlPts(:,1,3,2) = [  x_1+k	0           y_1-k	 1           ];
controlPts(:,2,3,2) = [  x_1+k	x_1+k       y_1-k	 1/sqrt(2)   ];
controlPts(:,3,3,2) = [  0      x_1+k       y_1-k    1           ];
controlPts(:,4,3,2) = [ -(x_1+k)   x_1+k   	y_1-k    1/sqrt(2)   ];
controlPts(:,5,3,2) = [ -(x_1+k)   0        y_1-k    1           ];
controlPts(:,6,3,2) = [ -(x_1+k)  -(x_1+k) 	y_1-k    1/sqrt(2)   ];
controlPts(:,7,3,2) = [  0     -(x_1+k)     y_1-k    1           ];
controlPts(:,8,3,2) = [  x_1+k  -(x_1+k)   	y_1-k    1/sqrt(2)   ];
controlPts(:,9,3,2) = [  x_1+k	0           y_1-k    1           ];

controlPts(:,1,4,2) = [  x_1	0        y_1	1           ];
controlPts(:,2,4,2) = [  x_1	x_1      y_1	1/sqrt(2)   ];
controlPts(:,3,4,2) = [  0      x_1      y_1    1           ];
controlPts(:,4,4,2) = [ -x_1    x_1      y_1    1/sqrt(2)   ];
controlPts(:,5,4,2) = [ -x_1    0        y_1    1           ];
controlPts(:,6,4,2) = [ -x_1   -x_1      y_1    1/sqrt(2)   ];
controlPts(:,7,4,2) = [  0     -x_1      y_1    1           ];
controlPts(:,8,4,2) = [  x_1   -x_1      y_1    1/sqrt(2)   ];
controlPts(:,9,4,2) = [  x_1	0        y_1    1           ];

controlPts(:,1,5,2) = [  x_1-k	0           y_1+k	 1           ];
controlPts(:,2,5,2) = [  x_1-k	x_1-k       y_1+k	 1/sqrt(2)   ];
controlPts(:,3,5,2) = [  0      x_1-k       y_1+k    1           ];
controlPts(:,4,5,2) = [ -(x_1-k)   x_1-k   	y_1+k    1/sqrt(2)   ];
controlPts(:,5,5,2) = [ -(x_1-k)   0        y_1+k    1           ];
controlPts(:,6,5,2) = [ -(x_1-k)  -(x_1-k) 	y_1+k    1/sqrt(2)   ];
controlPts(:,7,5,2) = [  0     -(x_1-k)     y_1+k    1           ];
controlPts(:,8,5,2) = [  x_1-k  -(x_1-k)   	y_1+k    1/sqrt(2)   ];
controlPts(:,9,5,2) = [  x_1-k	0           y_1+k    1           ];

controlPts(:,1,6,2) = [  R_4o	0        (h1+h2)/2	  1           ];
controlPts(:,2,6,2) = [  R_4o	R_4o     (h1+h2)/2	  1/sqrt(2)   ];
controlPts(:,3,6,2) = [  0      R_4o     (h1+h2)/2    1           ];
controlPts(:,4,6,2) = [ -R_4o    R_4o    (h1+h2)/2    1/sqrt(2)   ];
controlPts(:,5,6,2) = [ -R_4o    0       (h1+h2)/2    1           ];
controlPts(:,6,6,2) = [ -R_4o   -R_4o    (h1+h2)/2    1/sqrt(2)   ];
controlPts(:,7,6,2) = [  0     -R_4o     (h1+h2)/2    1           ];
controlPts(:,8,6,2) = [  R_4o   -R_4o    (h1+h2)/2    1/sqrt(2)   ];
controlPts(:,9,6,2) = [  R_4o	0        (h1+h2)/2    1           ];

controlPts(:,1,7,2) = [  x_2-k	0           y_2-k	 1           ];
controlPts(:,2,7,2) = [  x_2-k	x_2-k       y_2-k	 1/sqrt(2)   ];
controlPts(:,3,7,2) = [  0      x_2-k       y_2-k    1           ];
controlPts(:,4,7,2) = [ -(x_2-k)   x_2-k   	y_2-k    1/sqrt(2)   ];
controlPts(:,5,7,2) = [ -(x_2-k)   0        y_2-k    1           ];
controlPts(:,6,7,2) = [ -(x_2-k)  -(x_2-k) 	y_2-k    1/sqrt(2)   ];
controlPts(:,7,7,2) = [  0     -(x_2-k)     y_2-k    1           ];
controlPts(:,8,7,2) = [  x_2-k  -(x_2-k)   	y_2-k    1/sqrt(2)   ];
controlPts(:,9,7,2) = [  x_2-k	0           y_2-k    1           ];

controlPts(:,1,8,2) = [  x_2	0        y_2	1           ];
controlPts(:,2,8,2) = [  x_2	x_2      y_2	1/sqrt(2)   ];
controlPts(:,3,8,2) = [  0      x_2      y_2    1           ];
controlPts(:,4,8,2) = [ -x_2    x_2      y_2    1/sqrt(2)   ];
controlPts(:,5,8,2) = [ -x_2    0        y_2    1           ];
controlPts(:,6,8,2) = [ -x_2   -x_2      y_2    1/sqrt(2)   ];
controlPts(:,7,8,2) = [  0     -x_2      y_2    1           ];
controlPts(:,8,8,2) = [  x_2   -x_2      y_2    1/sqrt(2)   ];
controlPts(:,9,8,2) = [  x_2	0        y_2    1           ];

controlPts(:,1,9,2) = [  x_2+k	0           y_2+k	 1           ];
controlPts(:,2,9,2) = [  x_2+k	x_2+k       y_2+k	 1/sqrt(2)   ];
controlPts(:,3,9,2) = [  0      x_2+k       y_2+k    1           ];
controlPts(:,4,9,2) = [ -(x_2+k)   x_2+k   	y_2+k    1/sqrt(2)   ];
controlPts(:,5,9,2) = [ -(x_2+k)   0        y_2+k    1           ];
controlPts(:,6,9,2) = [ -(x_2+k)  -(x_2+k) 	y_2+k    1/sqrt(2)   ];
controlPts(:,7,9,2) = [  0     -(x_2+k)     y_2+k    1           ];
controlPts(:,8,9,2) = [  x_2+k  -(x_2+k)   	y_2+k    1/sqrt(2)   ];
controlPts(:,9,9,2) = [  x_2+k	0           y_2+k    1           ];

controlPts(:,1,10,2) = [  R_6o	0          h6	1           ];
controlPts(:,2,10,2) = [  R_6o	R_6o       h6	1/sqrt(2)   ];
controlPts(:,3,10,2) = [  0      R_6o      h6    1           ];
controlPts(:,4,10,2) = [ -R_6o    R_6o     h6    1/sqrt(2)   ];
controlPts(:,5,10,2) = [ -R_6o    0        h6    1           ];
controlPts(:,6,10,2) = [ -R_6o   -R_6o     h6    1/sqrt(2)   ];
controlPts(:,7,10,2) = [  0     -R_6o      h6    1           ];
controlPts(:,8,10,2) = [  R_6o   -R_6o     h6    1/sqrt(2)   ];
controlPts(:,9,10,2) = [  R_6o	0          h6    1           ];

controlPts(:,1,11,2) = [  R_7o	0          h7	1           ];
controlPts(:,2,11,2) = [  R_7o	R_7o       h7	1/sqrt(2)   ];
controlPts(:,3,11,2) = [  0      R_7o      h7    1           ];
controlPts(:,4,11,2) = [ -R_7o    R_7o     h7    1/sqrt(2)   ];
controlPts(:,5,11,2) = [ -R_7o    0        h7    1           ];
controlPts(:,6,11,2) = [ -R_7o   -R_7o     h7    1/sqrt(2)   ];
controlPts(:,7,11,2) = [  0     -R_7o      h7    1           ];
controlPts(:,8,11,2) = [  R_7o   -R_7o     h7    1/sqrt(2)   ];
controlPts(:,9,11,2) = [  R_7o	0          h7    1           ];

controlPts(:,1,12,2) = [  R_8o	0          h8	 1           ];
controlPts(:,2,12,2) = [  R_8o	R_8o       h8	 1/sqrt(2)   ];
controlPts(:,3,12,2) = [  0      R_8o      h8    1           ];
controlPts(:,4,12,2) = [ -R_8o    R_8o     h8    1/sqrt(2)   ];
controlPts(:,5,12,2) = [ -R_8o    0        h8    1           ];
controlPts(:,6,12,2) = [ -R_8o   -R_8o     h8    1/sqrt(2)   ];
controlPts(:,7,12,2) = [  0     -R_8o      h8    1           ];
controlPts(:,8,12,2) = [  R_8o   -R_8o     h8    1/sqrt(2)   ];
controlPts(:,9,12,2) = [  R_8o	0          h8    1           ];

controlPts(:,1,13,2) = [  R_8o	0          h9	 1           ];
controlPts(:,2,13,2) = [  R_8o	R_8o       h9	 1/sqrt(2)   ];
controlPts(:,3,13,2) = [  0      R_8o      h9    1           ];
controlPts(:,4,13,2) = [ -R_8o    R_8o     h9    1/sqrt(2)   ];
controlPts(:,5,13,2) = [ -R_8o    0        h9    1           ];
controlPts(:,6,13,2) = [ -R_8o   -R_8o     h9    1/sqrt(2)   ];
controlPts(:,7,13,2) = [  0     -R_8o      h9    1           ];
controlPts(:,8,13,2) = [  R_8o   -R_8o     h9    1/sqrt(2)   ];
controlPts(:,9,13,2) = [  R_8o	0          h9    1           ];

nurbs{1} = createNURBSobject(controlPts,{Xi, Eta, Zeta});
nurbs{1} = insertKnotsInNURBS(nurbs{1},{[] [0.5, 1.1, 1.9]/4 0.5});
nurbs{1} = elevateNURBSdegree(nurbs{1},[0 0 1]);



Xi = nurbs{1}.knots{1};
Eta = [0,0,0,0.5,0.5,1,1,1];
Zeta = [0,0,1,1];

controlPts = zeros(4,9,5,2);

controlPts(:,1,1,2) = [  R_0o	0       0	 1           ];
controlPts(:,2,1,2) = [  R_0o	R_0o    0	 1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0      R_0o    0    1           ];
controlPts(:,4,1,2) = [ -R_0o	R_0o    0    1/sqrt(2)   ];
controlPts(:,5,1,2) = [ -R_0o   0       0    1           ];
controlPts(:,6,1,2) = [ -R_0o  -R_0o    0    1/sqrt(2)   ];
controlPts(:,7,1,2) = [  0     -R_0o    0    1           ];
controlPts(:,8,1,2) = [  R_0o  -R_0o    0    1/sqrt(2)   ];
controlPts(:,9,1,2) = [  R_0o	0       0    1           ];


controlPts(:,1,2,2) = [  (R_0o+t2/2)	0           0	 1/sqrt(2)           ];
controlPts(:,2,2,2) = [  (R_0o+t2/2)	(R_0o+t2/2)  	0	 1/2   ];
controlPts(:,3,2,2) = [  0          (R_0o+t2/2) 	0    1/sqrt(2)           ];
controlPts(:,4,2,2) = [ -(R_0o+t2/2)  (R_0o+t2/2) 	0    1/2   ];
controlPts(:,5,2,2) = [ -(R_0o+t2/2)   0          0    1/sqrt(2)           ];
controlPts(:,6,2,2) = [ -(R_0o+t2/2)  -(R_0o+t2/2) 	0    1/2   ];
controlPts(:,7,2,2) = [  0          -(R_0o+t2/2)	0    1/sqrt(2)           ];
controlPts(:,8,2,2) = [  (R_0o+t2/2)  -(R_0o+t2/2)	0    1/2   ];
controlPts(:,9,2,2) = [  (R_0o+t2/2)	0           0    1/sqrt(2)           ];

controlPts(:,1,3,2) = [  (R_0o+t2/2)	0           t2/2	1           ];
controlPts(:,2,3,2) = [  (R_0o+t2/2)	(R_0o+t2/2)  	t2/2	1/sqrt(2)   ];
controlPts(:,3,3,2) = [  0          (R_0o+t2/2) 	t2/2    1           ];
controlPts(:,4,3,2) = [ -(R_0o+t2/2)  (R_0o+t2/2) 	t2/2    1/sqrt(2)   ];
controlPts(:,5,3,2) = [ -(R_0o+t2/2)   0          t2/2    1           ];
controlPts(:,6,3,2) = [ -(R_0o+t2/2)  -(R_0o+t2/2) 	t2/2    1/sqrt(2)   ];
controlPts(:,7,3,2) = [  0          -(R_0o+t2/2)	t2/2    1           ];
controlPts(:,8,3,2) = [  (R_0o+t2/2)  -(R_0o+t2/2)	t2/2    1/sqrt(2)   ];
controlPts(:,9,3,2) = [  (R_0o+t2/2)	0           t2/2    1           ];

controlPts(:,1,4,2) = [  (R_0o+t2/2)	0           t2      1/sqrt(2)           ];
controlPts(:,2,4,2) = [  (R_0o+t2/2)	(R_0o+t2/2)  	t2      1/2   ];
controlPts(:,3,4,2) = [  0          (R_0o+t2/2) 	t2      1/sqrt(2)           ];
controlPts(:,4,4,2) = [ -(R_0o+t2/2)  (R_0o+t2/2) 	t2      1/2   ];
controlPts(:,5,4,2) = [ -(R_0o+t2/2)   0          t2      1/sqrt(2)           ];
controlPts(:,6,4,2) = [ -(R_0o+t2/2)  -(R_0o+t2/2) 	t2      1/2   ];
controlPts(:,7,4,2) = [  0          -(R_0o+t2/2)	t2      1/sqrt(2)           ];
controlPts(:,8,4,2) = [  (R_0o+t2/2)  -(R_0o+t2/2)	t2      1/2   ];
controlPts(:,9,4,2) = [  (R_0o+t2/2)	0           t2      1/sqrt(2)           ];

controlPts(:,1,5,2) = [  R_0o	0        t2	  1           ];
controlPts(:,2,5,2) = [  R_0o	R_0o     t2	  1/sqrt(2)   ];
controlPts(:,3,5,2) = [  0      R_0o     t2    1           ];
controlPts(:,4,5,2) = [ -R_0o   R_0o     t2    1/sqrt(2)   ];
controlPts(:,5,5,2) = [ -R_0o   0        t2    1           ];
controlPts(:,6,5,2) = [ -R_0o  -R_0o     t2    1/sqrt(2)   ];
controlPts(:,7,5,2) = [  0     -R_0o     t2    1           ];
controlPts(:,8,5,2) = [  R_0o  -R_0o     t2    1/sqrt(2)   ];
controlPts(:,9,5,2) = [  R_0o	0        t2    1           ];

for i = 1:5
    controlPts(:,1,i,1) = [  R_0o	0        t2/2	  1           ];
    controlPts(:,2,i,1) = [  R_0o	R_0o     t2/2	  1/sqrt(2)   ];
    controlPts(:,3,i,1) = [  0      R_0o     t2/2    1           ];
    controlPts(:,4,i,1) = [ -R_0o   R_0o     t2/2    1/sqrt(2)   ];
    controlPts(:,5,i,1) = [ -R_0o   0        t2/2    1           ];
    controlPts(:,6,i,1) = [ -R_0o  -R_0o     t2/2    1/sqrt(2)   ];
    controlPts(:,7,i,1) = [  0     -R_0o     t2/2    1           ];
    controlPts(:,8,i,1) = [  R_0o  -R_0o     t2/2    1/sqrt(2)   ];
    controlPts(:,9,i,1) = [  R_0o	0        t2/2    1           ];
end
nurbs{2} = createNURBSobject(controlPts,{Xi, Eta, Zeta});
nurbs{2} = elevateNURBSdegree(nurbs{2},[0 0 1]);


controlPts(:,1,1,2) = [  R_8o	0       h9	 1           ];
controlPts(:,2,1,2) = [  R_8o	R_8o    h9	 1/sqrt(2)   ];
controlPts(:,3,1,2) = [  0      R_8o    h9   1           ];
controlPts(:,4,1,2) = [ -R_8o	R_8o    h9   1/sqrt(2)   ];
controlPts(:,5,1,2) = [ -R_8o   0       h9   1           ];
controlPts(:,6,1,2) = [ -R_8o  -R_8o    h9   1/sqrt(2)   ];
controlPts(:,7,1,2) = [  0     -R_8o    h9   1           ];
controlPts(:,8,1,2) = [  R_8o  -R_8o    h9   1/sqrt(2)   ];
controlPts(:,9,1,2) = [  R_8o	0       h9   1           ];


controlPts(:,1,2,2) = [  R_8o	0           h9+t1/2	 1/sqrt(2)           ];
controlPts(:,2,2,2) = [  R_8o	R_8o        h9+t1/2	 1/2   ];
controlPts(:,3,2,2) = [  0   	R_8o        h9+t1/2    1/sqrt(2)           ];
controlPts(:,4,2,2) = [ -R_8o  R_8o         h9+t1/2    1/2   ];
controlPts(:,5,2,2) = [ -R_8o   0           h9+t1/2    1/sqrt(2)           ];
controlPts(:,6,2,2) = [ -R_8o  -R_8o        h9+t1/2    1/2   ];
controlPts(:,7,2,2) = [  0   	-R_8o       h9+t1/2    1/sqrt(2)           ];
controlPts(:,8,2,2) = [  R_8o  -R_8o        h9+t1/2    1/2   ];
controlPts(:,9,2,2) = [  R_8o	0           h9+t1/2    1/sqrt(2)           ];

controlPts(:,1,3,2) = [  (R_8o+R_8i)/2	0               h9+t1/2	1           ];
controlPts(:,2,3,2) = [  (R_8o+R_8i)/2	(R_8o+R_8i)/2  	h9+t1/2	1/sqrt(2)   ];
controlPts(:,3,3,2) = [  0              (R_8o+R_8i)/2 	h9+t1/2    1           ];
controlPts(:,4,3,2) = [ -(R_8o+R_8i)/2  (R_8o+R_8i)/2 	h9+t1/2    1/sqrt(2)   ];
controlPts(:,5,3,2) = [ -(R_8o+R_8i)/2   0              h9+t1/2    1           ];
controlPts(:,6,3,2) = [ -(R_8o+R_8i)/2  -(R_8o+R_8i)/2 	h9+t1/2    1/sqrt(2)   ];
controlPts(:,7,3,2) = [  0              -(R_8o+R_8i)/2	h9+t1/2    1           ];
controlPts(:,8,3,2) = [  (R_8o+R_8i)/2  -(R_8o+R_8i)/2	h9+t1/2    1/sqrt(2)   ];
controlPts(:,9,3,2) = [  (R_8o+R_8i)/2	0               h9+t1/2    1           ];

controlPts(:,1,4,2) = [  R_8i	0   	h9+t1/2      1/sqrt(2)           ];
controlPts(:,2,4,2) = [  R_8i	R_8i  	h9+t1/2      1/2   ];
controlPts(:,3,4,2) = [  0   	R_8i 	h9+t1/2      1/sqrt(2)           ];
controlPts(:,4,4,2) = [ -R_8i   R_8i 	h9+t1/2      1/2   ];
controlPts(:,5,4,2) = [ -R_8i   0     	h9+t1/2      1/sqrt(2)           ];
controlPts(:,6,4,2) = [ -R_8i  -R_8i 	h9+t1/2      1/2   ];
controlPts(:,7,4,2) = [  0    	-R_8i	h9+t1/2      1/sqrt(2)           ];
controlPts(:,8,4,2) = [  R_8i  -R_8i	h9+t1/2      1/2   ];
controlPts(:,9,4,2) = [  R_8i	0       h9+t1/2      1/sqrt(2)           ];

controlPts(:,1,5,2) = [  R_8i	0        h9	  1           ];
controlPts(:,2,5,2) = [  R_8i	R_8i     h9	  1/sqrt(2)   ];
controlPts(:,3,5,2) = [  0      R_8i     h9    1           ];
controlPts(:,4,5,2) = [ -R_8i   R_8i     h9    1/sqrt(2)   ];
controlPts(:,5,5,2) = [ -R_8i   0        h9    1           ];
controlPts(:,6,5,2) = [ -R_8i  -R_8i     h9    1/sqrt(2)   ];
controlPts(:,7,5,2) = [  0     -R_8i     h9    1           ];
controlPts(:,8,5,2) = [  R_8i  -R_8i     h9    1/sqrt(2)   ];
controlPts(:,9,5,2) = [  R_8i	0        h9    1           ];

for i = 1:5
    controlPts(:,1,i,1) = [  (R_8o+R_8i)/2	0            	h9	  1           ];
    controlPts(:,2,i,1) = [  (R_8o+R_8i)/2	(R_8o+R_8i)/2 	h9	  1/sqrt(2)   ];
    controlPts(:,3,i,1) = [  0              (R_8o+R_8i)/2  	h9    1           ];
    controlPts(:,4,i,1) = [ -(R_8o+R_8i)/2   (R_8o+R_8i)/2	h9    1/sqrt(2)   ];
    controlPts(:,5,i,1) = [ -(R_8o+R_8i)/2   0              h9    1           ];
    controlPts(:,6,i,1) = [ -(R_8o+R_8i)/2  -(R_8o+R_8i)/2	h9    1/sqrt(2)   ];
    controlPts(:,7,i,1) = [  0              -(R_8o+R_8i)/2 	h9    1           ];
    controlPts(:,8,i,1) = [  (R_8o+R_8i)/2  -(R_8o+R_8i)/2 	h9    1/sqrt(2)   ];
    controlPts(:,9,i,1) = [  (R_8o+R_8i)/2	0               h9    1           ];
end
nurbs{3} = createNURBSobject(controlPts,{Xi, Eta, Zeta});
nurbs{3} = elevateNURBSdegree(nurbs{3},[0 0 1]);

nurbs = nurbs([2,1,3]);


