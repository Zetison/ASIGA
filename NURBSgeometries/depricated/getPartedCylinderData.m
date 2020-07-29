function solid = getPartedCylinderData(R_i, R_o, L, theta, phi)
error('Depricated use getCylinderData')


if phi >= pi
    error('The angle must be less than 180 degrees')
end
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];


x_1i = R_i*cos(theta);
z_1i = R_i*sin(theta);

x_2i = R_i*(cos(theta)-sin(theta)*tan(phi/2));
z_2i = R_i*(sin(theta)+cos(theta)*tan(phi/2));

x_3i = -R_i*cos(phi+theta);
z_3i = R_i*sin(phi+theta);

x_1o = R_o*cos(theta);
z_1o = R_o*sin(theta);

x_2o = R_o*(cos(theta)-sin(theta)*tan(phi/2));
z_2o = R_o*(sin(theta)+cos(theta)*tan(phi/2));

x_3o = -R_o*cos(phi+theta);
z_3o = R_o*sin(phi+theta);

weight = cos(phi/2);

controlPts = zeros(4,3,2,2);

controlPts(:,1,1,1) = [  x_1i	0    	z_1i    1           ];
controlPts(:,2,1,1) = [  x_2i 	0       z_2i    weight      ];
controlPts(:,3,1,1) = [  x_3i   0       z_3i    1           ];

controlPts(:,1,2,1) = [  x_1i	L    	z_1i    1           ];
controlPts(:,2,2,1) = [  x_2i 	L       z_2i    weight      ];
controlPts(:,3,2,1) = [  x_3i   L       z_3i    1           ];

controlPts(:,1,1,2) = [  x_1o	0    	z_1o    1           ];
controlPts(:,2,1,2) = [  x_2o 	0       z_2o    weight      ];
controlPts(:,3,1,2) = [  x_3o   0       z_3o    1           ];

controlPts(:,1,2,2) = [  x_1o	L    	z_1o    1           ];
controlPts(:,2,2,2) = [  x_2o 	L       z_2o    weight      ];
controlPts(:,3,2,2) = [  x_3o   L       z_3o    1           ];

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});

