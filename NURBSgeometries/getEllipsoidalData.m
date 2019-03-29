function solid = getEllipsoidalData(c_x,c_y,c_z,alignWithAxis, x_0, alpha, iXi, iEta)
if nargin < 4
    alignWithAxis = 'Xaxis';
end
if nargin < 5
    x_0 = [0, 0, 0];
end
if nargin < 6
    alpha = 0;
end
if nargin < 7
    iXi = [1,1,2,2,3,3]/4;
end
if nargin < 8
    iEta = [1,1]/2;
end

Xi = [0 0 0 iXi 1 1 1];
Eta = [0 0 0 iEta 1 1 1];

ctrlPtsXi = parmArc(Xi,2*pi);
ctrlPtsEta = parmArc(Eta,pi);

ctrlPtsEta(1,:) = -ctrlPtsEta(1,:); % flip parametrization (to go along x-axis)
controlPts = calcTensorRotCtrlPts(ctrlPtsXi,ctrlPtsEta);

R_x = rotationMatrix(alpha, 'Xaxis');
for i = 1:size(controlPts,2)
    for j = 1:size(controlPts,3)
        controlPts(1:3,i,j) = R_x*controlPts(1:3,i,j);
    end
end

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
controlPts(1,:,:) = c_x*controlPts(1,:,:) + x_0(1);
controlPts(2,:,:) = c_y*controlPts(2,:,:) + x_0(2);
controlPts(3,:,:) = c_z*controlPts(3,:,:) + x_0(3);


solid = createNURBSobject(controlPts,{Xi, Eta});


