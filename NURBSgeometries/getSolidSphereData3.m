function solid = getSolidSphereData3(R, alignWithAxis,www,p)
if nargin < 2
    alignWithAxis = 'Zaxis';
end
Xi = [0 0 0 0 1 1 1 1];
Eta = [0 0 0 0 1 1 1 1];
Zeta = [0 0 0 0 1 1 1 1];

controlPts = zeros(4,4,4,4);
L = R/sqrt(3);

% w = 1/sqrt(2);
theta = asin(1/sqrt(3));
w = 1/9*(3+2*sqrt(6));
% w = www(1);
w2 = w^2; % 0.770294776708619

% w2 = www(2);
% s = R/sqrt(6);
r = L*tan(theta);
s = r/sqrt(2); %= 1/2
L3 = L;
L4 = 2*(L+s);
% L4 = p;
L2 = sqrt(11-4*sqrt(6));
L5 = sqrt(29-6*sqrt(6));

x1 = 1.436364862287918;
x2 = 5.787524313268550;
% 
x1 = p(1);
x2 = p(2);
w2 = www;



controlPts(:,1,1,1) = [	-L	-L	-L	1]; % interpolatory point
controlPts(:,2,1,1) = [	 0  -L2	-L2	w];
controlPts(:,3,1,1) = [	 0  -L2	-L2	w];
controlPts(:,4,1,1) = [	 L	-L	-L	1]; % interpolatory point

controlPts(:,1,2,1) = [	-L2	 0  -L2	w];
controlPts(:,2,2,1) = [	 0       0      -L4	w2]; % ??
controlPts(:,3,2,1) = [	 0       0      -L4	w2]; % ??
controlPts(:,4,2,1) = [	 L2	 0      -L2	w];

controlPts(:,1,3,1) = [	-L2	 0  -L2	w];
controlPts(:,2,3,1) = [	 0   0  -L4	w2]; % ??
controlPts(:,3,3,1) = [	 0   0  -L4	w2]; % ??
controlPts(:,4,3,1) = [	 L2	 0  -L2	w];

controlPts(:,1,4,1) = [	-L	 L	-L	1]; % interpolatory point
controlPts(:,2,4,1) = [	 0   L2	-L2	w];
controlPts(:,3,4,1) = [	 0   L2	-L2	w];
controlPts(:,4,4,1) = [	 L	 L	-L	1]; % interpolatory point

controlPts(:,1,1,2) = [	-L-L3	-L-L3	 0  	1];
controlPts(:,2,1,2) = [	 0      -L-L3	 0  	w];
controlPts(:,3,1,2) = [	 0      -L-L3	 0  	w];
controlPts(:,4,1,2) = [	 L+L3	-L-L3	 0  	1];

controlPts(:,1,2,2) = [	-L-L3	 0       0  	1];
controlPts(:,2,2,2) = [	 0       0       0  	w];
controlPts(:,3,2,2) = [	 0       0       0  	w];
controlPts(:,4,2,2) = [	 L+L3	 0       0  	1];

controlPts(:,1,3,2) = [	-L-L3	 0       0  	1];
controlPts(:,2,3,2) = [	 0       0       0  	w];
controlPts(:,3,3,2) = [	 0       0       0  	w];
controlPts(:,4,3,2) = [	 L+L3	 0       0  	1];

controlPts(:,1,4,2) = [	-L-L3	 L+L3	 0  	1];
controlPts(:,2,4,2) = [	 0       L+L3	 0  	w];
controlPts(:,3,4,2) = [	 0       L+L3	 0  	w];
controlPts(:,4,4,2) = [	 L+L3	 L+L3	 0  	1];

controlPts(:,1,1,3) = [	-L-L3	-L-L3	 0  	1];
controlPts(:,2,1,3) = [	 0      -L-L3	 0  	w];
controlPts(:,3,1,3) = [	 0      -L-L3	 0  	w];
controlPts(:,4,1,3) = [	 L+L3	-L-L3	 0  	1];

controlPts(:,1,2,3) = [	-L-L3	 0       0  	1];
controlPts(:,2,2,3) = [	 0       0       0  	w];
controlPts(:,3,2,3) = [	 0       0       0  	w];
controlPts(:,4,2,3) = [	 L+L3	 0       0  	1];

controlPts(:,1,3,3) = [	-L-L3	 0       0  	1];
controlPts(:,2,3,3) = [	 0       0       0  	w];
controlPts(:,3,3,3) = [	 0       0       0  	w];
controlPts(:,4,3,3) = [	 L+L3	 0       0  	1];

controlPts(:,1,4,3) = [	-L-L3	 L+L3	 0  	1];
controlPts(:,2,4,3) = [	 0       L+L3	 0  	w];
controlPts(:,3,4,3) = [	 0       L+L3	 0  	w];
controlPts(:,4,4,3) = [	 L+L3	 L+L3	 0  	1];

controlPts(:,1,1,4) = [	-L	-L	L	1]; % interpolatory point
controlPts(:,2,1,4) = [	-L2  -L5	L5	w];
controlPts(:,3,1,4) = [	 L2  -L5	L5	w];
controlPts(:,4,1,4) = [	 L	-L	L	1]; % interpolatory point

controlPts(:,1,2,4) = [	-L5	 -L2  L5	w];
controlPts(:,2,2,4) = [	 -x1   -x1  x2	w2]; % ??
controlPts(:,3,2,4) = [	 x1   -x1  x2	w2]; % ??
controlPts(:,4,2,4) = [	L5	 -L2  L5	w];

controlPts(:,1,3,4) = [	-L5	 L2  L5	w];
controlPts(:,2,3,4) = [	 -x1   x1  x2	w2]; % ??
controlPts(:,3,3,4) = [	 x1   x1  x2	w2]; % ??
controlPts(:,4,3,4) = [	L5	 L2  L5	w];

controlPts(:,1,4,4) = [	-L	L	L	1]; % interpolatory point
controlPts(:,2,4,4) = [	-L2  L5	L5	w];
controlPts(:,3,4,4) = [	 L2  L5	L5	w];
controlPts(:,4,4,4) = [	 L	L	L	1]; % interpolatory point

switch alignWithAxis
    case 'Xaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = temp;
    case 'Yaxis' 
        temp = controlPts(1,:,:,:);
        controlPts(1,:,:,:) = controlPts(2,:,:,:);
        controlPts(2,:,:,:) = controlPts(3,:,:,:);
        controlPts(3,:,:,:) = temp;
    case 'Zaxis'
        % Nothing to be done
end



solid = createNURBSobject(controlPts,{Xi, Eta, Zeta});
