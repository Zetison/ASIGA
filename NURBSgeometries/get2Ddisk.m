function nurbs = get2Ddisk(R, translateX, Xi)
if nargin < 3
    Xi = [0 0 0 1 1 2 2 3 3 3]/3;
end
Eta = [0 0 1 1];
n = numel(Xi)-(2+1);
controlPts = zeros(4,n,2);
controlPts(4,:,:) = 1;
controlPts(2:4,:,2) = R*parmArc(Xi,2*pi);
controlPts(1,:,:) = translateX;

nurbs = createNURBSobject(controlPts,{Xi, Eta});
nurbs = rotateNURBS(nurbs,-pi/2,'Xaxis');