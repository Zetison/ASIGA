function solid = getWedgeData(R, t, theta, x, A)
if nargin < 3
    theta = pi/4;
end
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];
w = cos(theta/2);
controlPts = zeros(4,3,2,2);

controlPts(:,1,1,1) = [  0      0       0    1           ];
controlPts(:,2,1,1) = [  0      0       0    w           ];
controlPts(:,3,1,1) = [  0      0       0    1           ];

controlPts(:,1,2,1) = [  R      0       0    1           ];
controlPts(:,2,2,1) = [  R      R       0    w           ];
controlPts(:,3,2,1) = [  0   	R       0    1           ];

controlPts(:,1,1,2) = [  0      0       t    1           ];
controlPts(:,2,1,2) = [  0      0       t    w           ];
controlPts(:,3,1,2) = [  0      0       t    1           ];

controlPts(:,1,2,2) = [  R      0       t    1           ];
controlPts(:,2,2,2) = [  R      R       t    w           ];
controlPts(:,3,2,2) = [  0   	R       t    1           ];

for i = 1:size(controlPts,2)
    for j = 1:size(controlPts,3)
        for k = 1:size(controlPts,4)
            controlPts(1:3,i,j,k) = A*controlPts(1:3,i,j,k) + x;
        end
    end
end

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});