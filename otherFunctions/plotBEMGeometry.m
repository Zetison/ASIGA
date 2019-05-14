function pD = plotBEMGeometry(patches,plotGP)
pD.plotGP = plotGP;
if plotGP
    pD.plotPointsAsSpheres = 1;
    pD.pointsRadius = 6e-3;
    pD.noSpherePts = 20;
    pD.lineColor = 'blue';
    pD.lineStyle = '-';
    pD.patches = patches;
    if plotGP
        close all
        for patch = 1:numel(patches)
            plotNURBS(patches{patch}.nurbs,{'resolution',[100 100], 'elementBasedSamples',true,'samplingDistance',0.1});
        end
        axis equal
        axis off
        set(gca, 'Color', 'none');
    %     view(-100,20)
        view(10,20)
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