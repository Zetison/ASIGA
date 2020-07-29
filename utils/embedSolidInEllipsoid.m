function waterMesh = embedSolidInEllipsoid(nurbs,C,alignWithAxis,x_0, alpha, Xi, Eta)
if size(nurbs) > 1
    error('Not implemented')
end
if nargin < 5
    alpha = 0;
end
if nargin < 6 || any(isnan(Xi))
    Xi = nurbs{1}.knots{1};
end
if nargin < 7 || any(isnan(Eta))
    Eta = nurbs{1}.knots{2};
end

ellipsoid = getEllipsoidData('C',C,'alignWithAxis',alignWithAxis,'x_0',x_0, 'alpha', alpha, 'Xi', Xi, 'Eta', Eta);

Zeta = [0 0 1 1];
controlPts = nurbs{1}.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = ellipsoid{1}.coeffs;
          
waterMesh = createNURBSobject(controlPts,{Xi, Eta, Zeta});