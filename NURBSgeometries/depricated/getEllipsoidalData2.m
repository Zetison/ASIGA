function solid = getEllipsoidalData2(R,c_x,alignWithAxis, x_0, xi_trans, xi_trans2)
error('Depricated: use getEllipsoidData() instead')
% if nargin < 3
%     alignWithAxis = 'Xaxis';
% end
% if nargin < 5
%     xi_trans = 1/3;
% end
% Xi = [0, 0, 0, xi_trans, xi_trans, xi_trans2, xi_trans2, 1, 1, 1]/3;
% Eta = [0 0 0 1 1 2 2 2]/2;
% 
% controlPts = zeros(4,7,5);
% 
% % outer surface
% controlPts(:,1,1) = [ -c_x   0     0      1           ];
% controlPts(:,2,1) = [ -c_x   0     0      1/2   ];
% controlPts(:,3,1) = [ -c_x   0     0      1           ];
% controlPts(:,4,1) = [ -c_x   0     0      1/2   ];
% controlPts(:,5,1) = [ -c_x   0     0      1           ];
% controlPts(:,6,1) = [ -c_x   0     0      1/2   ];
% controlPts(:,7,1) = [ -c_x   0     0      1           ];
% 
% controlPts(:,1,2) = [ 	-c_x  R     0           1/sqrt(2)         ];
% controlPts(:,2,2) = [ 	-c_x  R     sqrt(3)*R  	1/2*1/sqrt(2)   ];
% controlPts(:,3,2) = [	-c_x -R/2   sqrt(3)/2*R	1/sqrt(2)         ];
% controlPts(:,4,2) = [   -c_x -2*R   0           1/2*1/sqrt(2)   ];
% controlPts(:,5,2) = [	-c_x -R/2  -sqrt(3)/2*R	1/sqrt(2)         ];
% controlPts(:,6,2) = [ 	-c_x  R    -sqrt(3)*R 	1/2*1/sqrt(2)   ];
% controlPts(:,7,2) = [  	-c_x  R     0           1/sqrt(2)         ];
% 
% controlPts(:,1,3) = [  	0  R     0           1   ];
% controlPts(:,2,3) = [	0  R     sqrt(3)*R  	1/2           ];
% controlPts(:,3,3) = [	0 -R/2   sqrt(3)/2*R	1   ];
% controlPts(:,4,3) = [	0 -2*R   0           1/2          ];
% controlPts(:,5,3) = [	0 -R/2  -sqrt(3)/2*R	1   ];
% controlPts(:,6,3) = [	0  R    -sqrt(3)*R 	1/2           ];
% controlPts(:,7,3) = [   0  R     0           1   ];
% 
% controlPts(:,1,4) = [	c_x  R     0           1/sqrt(2)         ];
% controlPts(:,2,4) = [ 	c_x  R     sqrt(3)*R  	1/2*1/sqrt(2)   ];
% controlPts(:,3,4) = [  	c_x -R/2   sqrt(3)/2*R	1/sqrt(2)       ];
% controlPts(:,4,4) = [	c_x -2*R   0           1/2*1/sqrt(2)   ];
% controlPts(:,5,4) = [ 	c_x -R/2  -sqrt(3)/2*R	1/sqrt(2)       ];
% controlPts(:,6,4) = [   c_x  R    -sqrt(3)*R 	1/2*1/sqrt(2)   ];
% controlPts(:,7,4) = [ 	c_x  R     0           1/sqrt(2)         ];
% 
% controlPts(:,1,5) = [  c_x  0     0       1           ];
% controlPts(:,2,5) = [  c_x  0     0     1/2   ];
% controlPts(:,3,5) = [  c_x  0     0      1           ];
% controlPts(:,4,5) = [  c_x  0     0       1/2   ];
% controlPts(:,5,5) = [  c_x  0     0       1           ];
% controlPts(:,6,5) = [  c_x  0     0       1/2   ];
% controlPts(:,7,5) = [  c_x  0     0        1           ];
% 
% switch alignWithAxis
%     case 'Xaxis' 
%         % Nothing to be done
%     case 'Yaxis' 
%         temp = controlPts(1,:,:,:);
%         controlPts(1,:,:,:) = controlPts(2,:,:,:);
%         controlPts(2,:,:,:) = controlPts(3,:,:,:);
%         controlPts(3,:,:,:) = temp;
%     case 'Zaxis'
%         temp = controlPts(1,:,:,:);
%         controlPts(1,:,:,:) = controlPts(2,:,:,:);
%         controlPts(2,:,:,:) = controlPts(3,:,:,:);
%         controlPts(3,:,:,:) = temp;
% end
% 
% controlPts(1,:,:,:) = controlPts(1,:,:,:) + x_0(1);
% controlPts(2,:,:,:) = controlPts(2,:,:,:) + x_0(2);
% controlPts(3,:,:,:) = controlPts(3,:,:,:) + x_0(3);
% 
% 
% solid = createNURBSobject(controlPts,{Xi, Eta});
