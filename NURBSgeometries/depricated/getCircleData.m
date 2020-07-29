function nurbs = getCircleData(R)
error('Depricated: getArcData')
if nargin == 0
    R = 1;
end
Xi = [0 0 0 1 1 2 2 3 3 4 4 4];
Xi = Xi/4;
controlPts = zeros(3,9);
controlPts(:,1) = [ 1   0   1           ];
controlPts(:,2) = [ 1   1   1/sqrt(2)   ];
controlPts(:,3) = [ 0   1   1           ];
controlPts(:,4) = [-1   1   1/sqrt(2)   ];
controlPts(:,5) = [-1   0   1           ];
controlPts(:,6) = [-1  -1   1/sqrt(2)   ];
controlPts(:,7) = [ 0  -1   1           ];
controlPts(:,8) = [ 1  -1   1/sqrt(2)   ];
controlPts(:,9) = [ 1   0   1           ];
controlPts(1:end-1,:) = controlPts(1:end-1,:)*R;
nurbs = createNURBSobject(controlPts,{Xi});