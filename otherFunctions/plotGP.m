function plotGP(pD,y,ptColor)
if plotPointsAsSpheres
    [Xs,Ys,Zs] = sphere(pD.noSpherePts);
    for ii = 1:size(y,1)
        surf(pD.pointsRadius*Xs+y(ii,1),pD.pointsRadius*Ys+y(ii,2),pD.pointsRadius*Zs+y(ii,3), 'FaceColor', ptColor, ...
                                                                                    'EdgeColor','none','LineStyle','none')
    end
else
    plot3(y(:,1),y(:,2),y(:,3),'*','color',ptColor)
end