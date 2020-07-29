function solid = getEllipsoidalData3(R,c_x,alignWithAxis, x_0)
error('Depricated: use getEllipsoidData() instead')
% if nargin < 3
%     alignWithAxis = 'Xaxis';
% end
% Xi = zeros(1,30);
% Xi(28:30) = 1;
% Xi(26:27) = 2/3;
% Xi(4:5) = 1/3;
% for i = 1:5
%     Xi(6+2*(i-1):7+2*(i-1)) = 1/3 + i/3/12;
%     Xi(24-2*(i-1):25-2*(i-1)) = 2/3 - i/3/12;
% end
% Eta = [0 0 0 1 1 2 2 2]/2;
% 
% t1 = 120*pi/180;
% p2 = t1/12;
% w1 = 0.5;
% w2 = cos(p2/2);
% w3 = cos(2*p2/2);
% cptsCircle = zeros(27,3);
% 
% p1 = t1-pi/2;
% cptsCircle(1:14,:) = [0,        -1,         1;
%                       NaN,      NaN,        w1;
%                       cos(p1),  sin(p1),    1;
%                       NaN,      NaN,        w2;
%                       cos(p1+p2),  sin(p1+p2),    1;
%                       NaN,      NaN,        w2;
%                       cos(p1+2*p2),  sin(p1+2*p2),    1;
%                       NaN,      NaN,        w2;
%                       cos(p1+3*p2),  sin(p1+3*p2),    1;
%                       NaN,      NaN,        w2;
%                       cos(p1+4*p2),  sin(p1+4*p2),    1;
%                       NaN,      NaN,        w2;
%                       cos(p1+5*p2),  sin(p1+5*p2),    1;
%                       NaN,      NaN,        w3];
% cptsCircle(15:end,:) = flipud(cptsCircle(1:13,:));
% cptsCircle(15:end,1) = -cptsCircle(15:end,1);
% for i = 1:13
%     y1 = cptsCircle(2*i-1,1);
%     z1 = cptsCircle(2*i-1,2);
%     y2 = cptsCircle(2*i+1,1);
%     z2 = cptsCircle(2*i+1,2);
%     detrm = y1*z2-y2*z1;
%     y3 = -z2^2*z1+(y1^2+z1^2)*z2-y2^2*z1;
%     z3 = -y2*y1^2+(y2^2+z2^2)*y1-y2*z1^2;
%     cptsCircle(2*i,1:2) = [y3,z3]/detrm;
% end
% 
% controlPts = zeros(4,27,5);
% for i = 1:5
%     controlPts(2:4,:,i) = cptsCircle.';
% end
% controlPts(2:3,:,[1,5]) = 0;
% controlPts(4,:,[2,4]) = controlPts(4,:,[2,4])/sqrt(2);
% 
% controlPts(1,:,1:2) = -1;
% controlPts(1,:,3) = 0;
% controlPts(1,:,4:5) = 1;
% 
% controlPts(1,:,:) = controlPts(1,:,:)*c_x;
% controlPts(2:3,:,:) = controlPts(2:3,:,:)*R;
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
