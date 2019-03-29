function solid = getRectangleData(L_x, L_y, centerCoordinate)


% If the centerCoordinate is not given, the prism is placed with its center
% at the origin

Xi = [0 0 1 1];
Eta = [0 0 1 1];

controlPts = zeros(3,2,2);


controlPts(:,1,1) = [	 -L_x/2	 -L_y/2	1];
controlPts(:,2,1) = [	 L_x/2	 -L_y/2	1];

controlPts(:,1,2) = [	 -L_x/2	 L_y/2	1];
controlPts(:,2,2) = [	 L_x/2	 L_y/2	1];


if nargin == 3
    for i = 1:2
        for j = 1:2
            controlPts(1:2,i,j) = controlPts(1:2,i,j) + centerCoordinate';
        end
    end
end

solid = createNURBSobject(controlPts,{Xi, Eta});
