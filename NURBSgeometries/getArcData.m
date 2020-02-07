function nurbs = getArcData(R)
if nargin == 0
    R = 1;
end
Xi = [0 0 0 1 1 1];
controlPts = zeros(3,3);
controlPts(:,1) = [ 1   0   1           ];
controlPts(:,2) = [ 1   1   1/sqrt(2)   ];
controlPts(:,3) = [ 0   1   1           ];
controlPts(1:end-1,:) = controlPts(1:end-1,:)*R;
nurbs = createNURBSobject(controlPts,{Xi});