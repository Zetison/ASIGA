function nurbs = getConeData(R_1, R_2, L, theta, phi, translateX)

if phi >= pi
    error('The angle must be less than 180 degrees')
end
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];


y1_1i = R_1;
z1_1i = 0;

y1_2i = R_1;
z1_2i = R_1*tan(phi/2);

y1_3i = R_1*cos(phi);
z1_3i = R_1*sin(phi);

y2_1i = R_2;
z2_1i = 0;

y2_2i = R_2;
z2_2i = R_2*tan(phi/2);

y2_3i = R_2*cos(phi);
z2_3i = R_2*sin(phi);

weight = cos(phi/2);

controlPts = zeros(4,3,2);

controlPts(:,1,1) = [  0, y1_1i     z1_1i    1           ];
controlPts(:,2,1) = [  0, y1_2i 	z1_2i    weight      ];
controlPts(:,3,1) = [  0, y1_3i     z1_3i    1           ];

controlPts(:,1,2) = [  L, y2_1i     z2_1i    1           ];
controlPts(:,2,2) = [  L, y2_2i     z2_2i    weight      ];
controlPts(:,3,2) = [  L, y2_3i  	z2_3i    1           ];

controlPts(1,:,:) = controlPts(1,:,:) + translateX;
R_x = rotationMatrix(theta, 'Xaxis');
for j = 1:2
    controlPts(1:3,:,j) = R_x*controlPts(1:3,:,j);
end

nurbs = createNURBSobject(controlPts,{Xi, Eta});

