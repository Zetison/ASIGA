function pD = plotBEMGeometry(nurbs,plotGP,npts,plotPointsAsSpheres,alphaValue)
if nargin < 5
    alphaValue = 1;
end
pD.plotGP = plotGP;
if plotGP
    pD.plotPointsAsSpheres = plotPointsAsSpheres;
    pD.pointsRadius = 6e-3;
    pD.noSpherePts = 20;
    pD.lineColor = 'blue';
    pD.lineStyle = '-';
    pD.nurbs = nurbs;
    if plotGP
        close all
%         for patch = 1:17
%         for patch = [189,191,220:228]
%         plotNURBS(nurbs,{'resolution',[npts npts], 'elementBasedSamples',true,'samplingDistance',0.1,'alphaValue',alphaValue});
        plotNURBS(nurbs,{'resolution',[npts npts],'alphaValue',alphaValue});
        pD.h = gca;
        axis equal
        axis off
        set(gca, 'Color', 'none');
%         view(-100,20) % BeTSSi side part
        view(10,20) % getTask_articleBEM_S1_quadVis
%         view(0,-90)
%         view(-90,0) % BeTSSi rear part
        drawnow
        hold on
%         if false
%             cp = zeros(size(cp_p,1),3);
%             for j = 1:size(cp_p,1)
%                 patch = patchIdx(j);
%                 cp(j,:) = evaluateNURBS(patches{patch}.nurbs, cp_p(j,:));
%             end
%             plotGP(pD,cp,'blue')
%         end
        ax = gca;               % get the current axis
        ax.Clipping = 'off';    % turn clipping off
        camlight
    end
end