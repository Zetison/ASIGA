function solid = getPartOfSphere(R, alignWithAxis,www,p)
if nargin < 2
    alignWithAxis = 'Zaxis';
end
Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

controlPts = zeros(4,3,3);
L = R/sqrt(3);
L2 = R/sqrt(2);

% w = 1/sqrt(2);
theta = asin(1/sqrt(3));
w = cos(pi/4/2);
w3 = cos(acos((1/3)*sqrt(3)*sqrt(2))/2);
% w = 1/sqrt(2);
w2 = www(1);
% w2 = w^2;
% s = R/sqrt(6);
r = L*tan(theta);
s = r/sqrt(2); %= 1/2
L3 = L;
L4 = (sqrt(3)-sqrt(2))*R;
% L4 = p;

controlPts(:,1,1) = [	0	0	 R	1]; % ok
controlPts(:,2,1) = [	R*tan(pi/8) 0  	 R	w]; % ok
controlPts(:,3,1) = [	 L2	0	 L2	 1]; % ok

controlPts(:,1,2) = [	0	 R*tan(pi/8)       R	w]; % ok
controlPts(:,2,2) = [	 p(1)      p(1)       R	w2];
controlPts(:,3,2) = [	 L2	 L4       L2	w3];

controlPts(:,1,3) = [	0	 L2	 L2	1]; % ok
controlPts(:,2,3) = [	 L4   L2	 L2	w3];
controlPts(:,3,3) = [	 L	 L	 L	1]; % ok

switch alignWithAxis
    case 'Xaxis' 
        temp = controlPts(1,:,:);
        controlPts(1,:,:) = controlPts(3,:,:);
        controlPts(3,:,:) = controlPts(2,:,:);
        controlPts(2,:,:) = temp;
    case 'Yaxis' 
        temp = controlPts(1,:,:);
        controlPts(1,:,:) = controlPts(2,:,:);
        controlPts(2,:,:) = controlPts(3,:,:);
        controlPts(3,:,:) = temp;
    case 'Zaxis'
        % Nothing to be done
end



solid = createNURBSobject(controlPts,{Xi, Eta});
