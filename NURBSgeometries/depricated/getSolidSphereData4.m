function solid = getSolidSphereData4(R, alignWithAxis,www,p)
error('Use getEllipsoidData() instead')
if nargin < 2
    alignWithAxis = 'Zaxis';
end
Xi = [0 0 0 0 0 1 1 1 1 1];
Eta = [0 0 0 0 0 1 1 1 1 1];
Zeta = [0 0 0 0 0 1 1 1 1 1];

controlPts = zeros(4,5,5,5);
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
% x1 = p(1);
% x2 = p(2);
% w2 = www;



controlPts(:,1,1,1) = [	-L	-L	-L	1]; % interpolatory point
controlPts(:,2,1,1) = [	 0  -L2	-L2	w];
controlPts(:,3,1,1) = [	 0  -L2	-L2	w];
controlPts(:,3,1,1) = [	 0  -L2	-L2	w];
controlPts(:,4,1,1) = [	 L	-L	-L	1]; % interpolatory point

controlPts(:,1,2,1) = [	-L2	 0  -L2	w];
controlPts(:,2,2,1) = [	 0       0      -L4	w2]; % ??
controlPts(:,3,2,1) = [	 0       0      -L4	w2]; % ??
controlPts(:,3,2,1) = [	 0       0      -L4	w2]; % ??
controlPts(:,4,2,1) = [	 L2	 0      -L2	w];

controlPts(:,1,3,1) = [	-L2	 0  -L2	w];
controlPts(:,2,3,1) = [	 0   0  -L4	w2]; % ??
controlPts(:,3,3,1) = [	 0   0  -L4	w2]; % ??
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

q1 = 1.589186225978911;
q2 = 3.535533905932738;
q3 = 3.781938860333538;

W1 = 0.908248290463863;
W2 = 0.877664387285151;
% w1 = 0.824914957130530;
% w2 = 0.797137179352752;
% w3 = 0.770294776708619;
w1 = www(1);
w2 = www(2);
w3 = www(3);

x1 = p(1);
x2 = p(2);
x3 = p(3);
x4 = p(4);
x5 = p(5);

controlPts(:,1,1,5) = [	-L	-L	L	1]; % interpolatory point
controlPts(:,2,1,5) = [	-q1  -q2	q2	W1];
controlPts(:,3,1,5) = [	 0   -q3	q3	W2];
controlPts(:,4,1,5) = [	 q1  -q2	q2	W1];
controlPts(:,5,1,5) = [	 L	-L	L	1]; % interpolatory point

controlPts(:,1,2,5) = [	-q2	 -L2  q2	W1];
controlPts(:,2,2,5) = [	 -x1  -x1  x2	w1]; % ??
controlPts(:,3,2,5) = [	 0    -x3  x4	w2]; % ??
controlPts(:,4,2,5) = [	 x1   -x1  x2	w1]; % ??
controlPts(:,5,2,5) = [	q2	 -L2  q2	W1];

controlPts(:,1,3,5) = [	-q3	 0  q3	W2];
controlPts(:,2,3,5) = [	 -x3   0  x4	w1]; % ??
controlPts(:,3,3,5) = [	 0    0   x5	w3]; % ??
controlPts(:,4,3,5) = [	 x3   0   x4	w1]; % ??
controlPts(:,5,3,5) = [	q3	 0  q3	W2];

controlPts(:,1,4,5) = [	-q2	 L2  q2	W1];
controlPts(:,2,4,5) = [	-x1   x1  x2	w1]; % ??
controlPts(:,3,4,5) = [	 0    x3  x4	w2]; % ??
controlPts(:,4,4,5) = [	 x1   x1  x2	w1]; % ??
controlPts(:,5,4,5) = [	q2	 L2  q2	W1];

controlPts(:,1,5,5) = [	-L	L	L	1]; % interpolatory point
controlPts(:,2,5,5) = [	-q1  q2	q2	W1];
controlPts(:,3,5,5) = [	 0   q3	q3	W2];
controlPts(:,4,5,5) = [	 q1  q2	q2	W1];
controlPts(:,5,5,5) = [	 L	L	L	1]; % interpolatory point

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
