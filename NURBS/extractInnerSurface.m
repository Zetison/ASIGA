function [nurbsSurf, varColSurf] = extractInnerSurface(nurbs, varCol)

% 
% if ~iscell(varCol)
%     varCol = {};
% end




% function [nurbsSurf, varColSurf] = extractInnerSurfacePatch(nurbs, varCol)
Xi = nurbs.knots{1};
Eta = nurbs.knots{2};

controlPts = nurbs.coeffs(:,:,:,1);

nurbsSurf = createNURBSobject(controlPts,{Xi, Eta});

if nargin == 2
    varColSurf = varCol;

    noSurfPts = numel(controlPts(1,:,:));
    varColSurf.weights = varCol.weights(1:noSurfPts);
    varColSurf.controlPts = varCol.controlPts(1:noSurfPts,:);
    varColSurf.nurbs = nurbsSurf;
    varColSurf.noCtrlPts = noSurfPts;
    varColSurf.noDofs = noSurfPts;

    varColSurf = generateIGA2DMesh_new(varColSurf);
    varColSurf = findDofsToRemove(varColSurf);
end