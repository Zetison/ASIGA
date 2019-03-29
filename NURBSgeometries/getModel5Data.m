function nurbs = getModel5Data(R, eta1, eta2, L, l, alignWithAxis)


Xi = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
Eta = [0, 0, eta1, eta2, 1, 1, eta1+1, eta2+1, 2, 2]/2;

controlPts1 = zeros(4,9,4);

% outer surface
controlPts1(:,1,1) = [ -L/2   0     0      1           ];
controlPts1(:,2,1) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts1(:,3,1) = [ -L/2   0     0      1           ];
controlPts1(:,4,1) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts1(:,5,1) = [ -L/2   0     0      1           ];
controlPts1(:,6,1) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts1(:,7,1) = [ -L/2   0     0      1           ];
controlPts1(:,8,1) = [ -L/2   0     0      1/sqrt(2)   ];
controlPts1(:,9,1) = [ -L/2   0     0      1           ];

controlPts1(:,1,2) = [ -L/2     R   0    	1           ];
controlPts1(:,2,2) = [ -L/2     R   R   	1/sqrt(2)   ];
controlPts1(:,3,2) = [	-L/2     0     R  	1           ];
controlPts1(:,4,2) = [	-L/2    -R   R  	1/sqrt(2)   ];
controlPts1(:,5,2) = [	-L/2    -R   0    	1           ];
controlPts1(:,6,2) = [	-L/2    -R  -R   	1/sqrt(2)   ];
controlPts1(:,7,2) = [	-L/2     0    -R   	1           ];
controlPts1(:,8,2) = [	-L/2     R  -R   	1/sqrt(2)   ];
controlPts1(:,9,2) = [	-L/2     R   0    	1           ];

controlPts1(:,1,3) = [  L/2     R   0    	1           ];
controlPts1(:,2,3) = [  L/2     R   R   	1/sqrt(2)   ];
controlPts1(:,3,3) = [	L/2     0     R  	1           ];
controlPts1(:,4,3) = [	L/2    -R   R  	1/sqrt(2)   ];
controlPts1(:,5,3) = [	L/2    -R   0    	1           ];
controlPts1(:,6,3) = [	L/2    -R  -R   	1/sqrt(2)   ];
controlPts1(:,7,3) = [	L/2     0    -R   	1           ];
controlPts1(:,8,3) = [	L/2     R  -R   	1/sqrt(2)   ];
controlPts1(:,9,3) = [	L/2     R   0    	1           ];

controlPts1(:,1,4) = [  L/2   0     0      1           ];
controlPts1(:,2,4) = [  L/2   0     0      1/sqrt(2)   ];
controlPts1(:,3,4) = [  L/2   0     0      1           ];
controlPts1(:,4,4) = [  L/2   0     0      1/sqrt(2)   ];
controlPts1(:,5,4) = [  L/2   0     0      1           ];
controlPts1(:,6,4) = [  L/2   0     0      1/sqrt(2)   ];
controlPts1(:,7,4) = [  L/2   0     0      1           ];
controlPts1(:,8,4) = [  L/2   0     0      1/sqrt(2)   ];
controlPts1(:,9,4) = [  L/2   0     0      1           ];

switch alignWithAxis
    case 'Xaxis' 
        % Nothing to be done
    case 'Yaxis' 
        temp = controlPts1(1,:,:,:);
        controlPts1(1,:,:,:) = controlPts1(2,:,:,:);
        controlPts1(2,:,:,:) = controlPts1(3,:,:,:);
        controlPts1(3,:,:,:) = temp;
    case 'Zaxis'
        temp = controlPts1(1,:,:,:);
        controlPts1(1,:,:,:) = controlPts1(2,:,:,:);
        controlPts1(2,:,:,:) = controlPts1(3,:,:,:);
        controlPts1(3,:,:,:) = temp;
end

% keyboard
switch alignWithAxis
    case 'Xaxis' 
        controlPts2 = controlPts1;
        controlPts2(3,:,:) = controlPts1(3,:,:) + l/2;
        controlPts1(3,:,:) = controlPts1(3,:,:) - l/2;
        controlPts = zeros(size(controlPts1).*[1, 1, 2]);
        controlPts(:,:,1:size(controlPts1,3)) = controlPts1;
        controlPts(:,:,size(controlPts1,3)+1:end) = controlPts2;
    case 'Zaxis'
        controlPts2 = controlPts1;
        controlPts2(1,:,:) = controlPts1(1,:,:) + l/2;
        controlPts1(1,:,:) = controlPts1(1,:,:) - l/2;
        controlPts = zeros(size(controlPts1).*[1, 1, 2]);
        controlPts(:,:,1:size(controlPts1,3)) = controlPts1;
        controlPts(:,:,size(controlPts1,3)+1:end) = controlPts2;
end

nurbs = createNURBSobject(controlPts,{Xi, Eta});
nurbs = elevateNURBSdegree(nurbs,[0 1]);
