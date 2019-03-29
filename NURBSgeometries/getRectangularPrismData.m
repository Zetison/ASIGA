function solid = getRectangularPrismData(L_x, L_y, L_z, centerCoordinate)


% If the centerCoordinate is not given, the prism is placed with its center
% at the origin

Xi = [0 0 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

controlPts = zeros(4,2,2,2);


controlPts(:,1,1,1) = [	 -L_x/2	 -L_y/2	-L_z/2	1];
controlPts(:,2,1,1) = [	 L_x/2	 -L_y/2	-L_z/2	1];

controlPts(:,1,2,1) = [	 -L_x/2	 L_y/2	-L_z/2	1];
controlPts(:,2,2,1) = [	 L_x/2	 L_y/2	-L_z/2	1];

controlPts(:,1,1,2) = [	 -L_x/2	 -L_y/2	L_z/2	1];
controlPts(:,2,1,2) = [	 L_x/2	 -L_y/2	L_z/2	1];

controlPts(:,1,2,2) = [	 -L_x/2	 L_y/2	L_z/2	1];
controlPts(:,2,2,2) = [	 L_x/2	 L_y/2	L_z/2	1];

if nargin == 4
    for i = 1:2
        for j = 1:2
            for k = 1:2
                controlPts(1:3,i,j,k) = controlPts(1:3,i,j,k) + centerCoordinate;
            end
        end
    end
end

solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
