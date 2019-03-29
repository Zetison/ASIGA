function nurbs = getTubeData(R, L, x0,alignWithAxis)
if nargin < 4
    alignWithAxis = 'Xaxis';
end

Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0,0,1,1];

controlPts = zeros(4,9,2);

controlPts(:,1,1) = [  0     R   0    	1           ];
controlPts(:,2,1) = [  0     R   R   	1/sqrt(2)   ];
controlPts(:,3,1) = [  0     0   R  	1           ];
controlPts(:,4,1) = [  0    -R   R  	1/sqrt(2)   ];
controlPts(:,5,1) = [  0    -R   0    	1           ];
controlPts(:,6,1) = [	0   -R  -R   	1/sqrt(2)   ];
controlPts(:,7,1) = [	0    0  -R   	1           ];
controlPts(:,8,1) = [	0    R  -R   	1/sqrt(2)   ];
controlPts(:,9,1) = [	0    R   0    	1           ];

controlPts(:,1,2) = [  L     R   0    	1           ];
controlPts(:,2,2) = [  L     R   R   	1/sqrt(2)   ];
controlPts(:,3,2) = [  L     0   R  	1           ];
controlPts(:,4,2) = [  L    -R   R  	1/sqrt(2)   ];
controlPts(:,5,2) = [  L    -R   0    	1           ];
controlPts(:,6,2) = [	L   -R  -R   	1/sqrt(2)   ];
controlPts(:,7,2) = [	L    0  -R   	1           ];
controlPts(:,8,2) = [	L    R  -R   	1/sqrt(2)   ];
controlPts(:,9,2) = [	L    R   0    	1           ];

switch alignWithAxis
    case 'Xaxis' 
        % Nothing to be done
    case 'Yaxis' 
        temp = controlPts(1,:,:);
        controlPts(1,:,:) = controlPts(3,:,:);
        controlPts(3,:,:) = controlPts(2,:,:);
        controlPts(2,:,:) = temp;
    case 'Zaxis'
        temp = controlPts(1,:,:);
        controlPts(1,:,:) = controlPts(2,:,:);
        controlPts(2,:,:) = controlPts(3,:,:);
        controlPts(3,:,:) = temp;
end
nurbs = createNURBSobject(controlPts,{Xi, Eta});
nurbs = translateNURBS(nurbs,x0);
