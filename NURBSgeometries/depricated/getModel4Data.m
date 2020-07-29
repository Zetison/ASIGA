function solid = getModel4Data(R, t, eta1, eta2)

error('Depricated use getBeTSSiM4Data')

if nargin < 3
    eta1 = R/(2*R+t);
    eta2 = (R+t)/(2*R+t);
end

if 0
    w5 = 0.806806632310317;
    c_1 = 0.968134681494708;
    c_2 = 0.285927513325901;
    Xi = [0 0 0 1 1 2 2 3 3 3]/3;
    Eta = [0 0 0 eta1 eta1 eta2 eta2 1 1 1];
    % 
    v_center = R/sqrt(3)*[1 1 1];
    v_bot = R/sqrt(2)*[1 1 0];
    v_west = R/sqrt(2)*[1 0 1];
    v_east = R/sqrt(2)*[0 1 1];
    v_Bot = c_1*v_bot + c_2*v_center;
    v_West = c_1*v_west + c_2*v_center;
    v_East = c_1*v_east + c_2*v_center;
    
    R2 = sqrt(R^2-2*t^2);
    if abs((R2/t)^2-1) < 10*eps
        R3 = t;
        R2 = t;
        theta2 = 0;
    else
        R3 = (-t-(t/R2)^2*t+R2+t^2/R2)/(1-(t/R2)^2);
        theta2 = pi/2-2*atan(t/R2);
    end
%     t2 = (R2-R)*R2/(sqrt(2)*t) + sqrt(2)*t;
    s = t + (R^2-2*t^2-R*sqrt(R^2-2*t^2))/(2*t);
    
    theta = asin(t/R);
    w = cos(theta/2);
    theta3 = asin(sqrt(2)*t/R);
    w2 = cos(theta2/2);
    w3 = cos(theta3/2);
%     w3*w
    controlPts = zeros(4,7,7);
    % outer surface
    controlPts(:,1,1) = [ 0   0     0      1           ];
    controlPts(:,2,1) = [ 0   0     0      1/sqrt(2)   ];
    controlPts(:,3,1) = [ 0   0     0      1           ];
    controlPts(:,4,1) = [ 0   0     0      1/sqrt(2)   ];
    controlPts(:,5,1) = [ 0   0     0      1           ];
    controlPts(:,6,1) = [ 0   0     0      1/sqrt(2)   ];
    controlPts(:,7,1) = [ 0   0     0      1           ];

    controlPts(:,1,2) = [  	R/2  0     0           1   ];
    controlPts(:,2,2) = [	R/2  R/2     0  	1/sqrt(2)           ];
    controlPts(:,3,2) = [	0  R/2     0	1   ];
    controlPts(:,4,2) = [	0  R/2    R/2           1/sqrt(2)          ];
    controlPts(:,5,2) = [	0  0    R/2	1   ];
    controlPts(:,6,2) = [	R/2  0    R/2 	1/sqrt(2)           ];
    controlPts(:,7,2) = [   R/2  0    0           1   ];

    controlPts(:,1,3) = [  	R  0     0           1   ];
    controlPts(:,2,3) = [	R  R     0  	1/sqrt(2)           ];
    controlPts(:,3,3) = [	0  R     0	1   ];
    controlPts(:,4,3) = [	0  R    R           1/sqrt(2)          ];
    controlPts(:,5,3) = [	0  0    R	1   ];
    controlPts(:,6,3) = [	R  0    R 	1/sqrt(2)           ];
    controlPts(:,7,3) = [   R  0    0           1   ];

    controlPts(:,1,4) = [  	R  s     s           w3   ];
    controlPts(:,2,4) = [	v_Bot  	w5           ];
    controlPts(:,3,4) = [	s  R     s      w3   ];
    controlPts(:,4,4) = [	v_East           w5         ];
    controlPts(:,5,4) = [	s  s    R     w3   ];
    controlPts(:,6,4) = [	v_West      w5           ];
    controlPts(:,7,4) = [   R  s    s           w3   ];

    controlPts(:,1,5) = [  	R2  t     t           1   ];
    controlPts(:,2,5) = [	R3  R3     t  	w2           ];
    controlPts(:,3,5) = [	t  R2     t	1   ];
    controlPts(:,4,5) = [	t  R3    R3           w2         ];
    controlPts(:,5,5) = [	t  t    R2	1   ];
    controlPts(:,6,5) = [	R3  t    R3 	w2           ];
    controlPts(:,7,5) = [   R2  t    t           1   ];

    controlPts(:,1,6) = [  	(R2+t)/2  t     t           1   ];
    controlPts(:,2,6) = [	(R3+t)/2  (R3+t)/2     t  	w2           ];
    controlPts(:,3,6) = [	t  (R2+t)/2     t	1   ];
    controlPts(:,4,6) = [	t  (R3+t)/2    (R3+t)/2           w2         ];
    controlPts(:,5,6) = [	t  t    (R2+t)/2	1   ];
    controlPts(:,6,6) = [	(R3+t)/2  t    (R3+t)/2 	w2           ];
    controlPts(:,7,6) = [   (R2+t)/2  t    t           1   ];

    controlPts(:,1,7) = [  t  t     t       1           ];
    controlPts(:,2,7) = [  t  t     t     1/sqrt(2)   ];
    controlPts(:,3,7) = [  t  t     t      1           ];
    controlPts(:,4,7) = [  t  t     t       1/sqrt(2)   ];
    controlPts(:,5,7) = [  t  t     t       1           ];
    controlPts(:,6,7) = [  t  t     t       1/sqrt(2)   ];
    controlPts(:,7,7) = [  t  t     t        1           ];
%     keyboard
else
    Xi = [0 0 0 1 1 2 2 3 3 3]/3;
    Eta = [0 0 eta1 eta2 1 1];
    
    controlPts = zeros(4,7,4);
    % outer surface
    controlPts(:,1,1) = [  -t/2  -t/2     -t/2  	1           ];
    controlPts(:,2,1) = [  -t/2  -t/2     -t/2      1/sqrt(2)   ];
    controlPts(:,3,1) = [  -t/2  -t/2     -t/2   	1           ];
    controlPts(:,4,1) = [  -t/2  -t/2     -t/2   	1/sqrt(2)   ];
    controlPts(:,5,1) = [  -t/2  -t/2     -t/2    	1           ];
    controlPts(:,6,1) = [  -t/2  -t/2     -t/2   	1/sqrt(2)   ];
    controlPts(:,7,1) = [  -t/2  -t/2     -t/2    	1           ];
    
    controlPts(:,1,2) = [  	R       -t/2 	-t/2  	1           ];
    controlPts(:,2,2) = [	R       R       -t/2  	1/sqrt(2)	];
    controlPts(:,3,2) = [	-t/2	R       -t/2 	1           ];
    controlPts(:,4,2) = [	-t/2	R       R     	1/sqrt(2) 	];
    controlPts(:,5,2) = [	-t/2	-t/2   	R       1           ];
    controlPts(:,6,2) = [	R       -t/2  	R       1/sqrt(2)  	];
    controlPts(:,7,2) = [   R      -t/2     -t/2   	1           ];

    controlPts(:,1,3) = [  	R   t/2 	t/2    	1           ];
    controlPts(:,2,3) = [	R   R       t/2  	1/sqrt(2)	];
    controlPts(:,3,3) = [	t/2	R       t/2     1           ];
    controlPts(:,4,3) = [	t/2	R       R     	1/sqrt(2) 	];
    controlPts(:,5,3) = [	t/2	t/2     R       1           ];
    controlPts(:,6,3) = [	R   t/2     R       1/sqrt(2)  	];
    controlPts(:,7,3) = [   R   t/2     t/2   	1           ];

    controlPts(:,1,4) = [  t/2  t/2     t/2  	1           ];
    controlPts(:,2,4) = [  t/2  t/2     t/2     1/sqrt(2)   ];
    controlPts(:,3,4) = [  t/2  t/2     t/2   	1           ];
    controlPts(:,4,4) = [  t/2  t/2     t/2   	1/sqrt(2)   ];
    controlPts(:,5,4) = [  t/2  t/2     t/2    	1           ];
    controlPts(:,6,4) = [  t/2  t/2     t/2   	1/sqrt(2)   ];
    controlPts(:,7,4) = [  t/2  t/2     t/2    	1           ];
end


solid = createNURBSobject(controlPts,{Xi, Eta});


