function waterMesh = embedSolidInEllipsoid(nurbs,a,b,c,alignWithAxis,x_0, alpha, iXi, iEta)

if nargin < 7
    alpha = 0;
end
Xi = nurbs.knots{1};
Eta = nurbs.knots{2};
if nargin < 8 || any(isnan(iXi))
    p_xi = nurbs.degree(1);
    iXi = Xi(p_xi+2:end-p_xi-1);
end
if nargin < 9 || any(isnan(iEta))
    p_eta = nurbs.degree(2);
    iEta = Eta(p_eta+2:end-p_eta-1);
end

ellipsoid = getEllipsoidalData(a,b,c,alignWithAxis,x_0, alpha, iXi, iEta);

Zeta = [0 0 1 1];
controlPts = nurbs.coeffs(:,:,:,end);
controlPts(:,:,:,end+1) = ellipsoid.coeffs;
          
waterMesh = createNURBSobject(controlPts,{Xi, Eta, Zeta});