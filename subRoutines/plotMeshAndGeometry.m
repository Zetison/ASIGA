
if strcmp(coreMethod, 'linear_FEM')
    prePlot.resolution = [0,0,0];
end
if prePlot.plot3Dgeometry
    figure('Color','white','name',['3D plot of geometry with mesh ' num2str(M)])
    if strcmp(method, 'IENSG')
        c_z = varCol{1}.c_z;
        c_xy = varCol{1}.c_xy;
        alignWithAxis = varCol{1}.alignWithAxis;
        x_0 = varCol{1}.x_0;
        ellipsoid = getEllipsoidData('C',[c_z,c_xy,c_xy],'alignWithAxis', alignWithAxis, 'x_0', x_0);
        alphaValue = 0.6;
        if prePlot.alphaValue == 1
            prePlot.alphaValue = 0.8;
        end
        plotNURBS(ellipsoid,'resolution',[20 40],'alphaValue',0.6,'color','blue');
    end
    for j = 1:numel(varCol)
        nurbs = varCol{j}.nurbs;
        prePlot.displayName = ['Domain ' num2str(j)];
        plotNURBS(nurbs, prePlot);
    end
    axis equal
    axis(prePlot.axis)
    if ~strcmp(prePlot.axis,'off')
        xlabel(prePlot.xlabel)
        ylabel(prePlot.ylabel)
        zlabel(prePlot.zlabel)
    end
    set(gca, 'Color', 'none');
    if ~isempty(prePlot.title)
        title(prePlot.title)
    end
    drawnow
    camproj('perspective')
%     camproj('orthographic')
    
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    if prePlot.useCamlight
        camlight
        material dull
    end
    if ~isempty(prePlot.export_fig_name3D)
        figureFullScreen(gcf)
        export_fig(prePlot.export_fig_name3D, '-png', '-transparent', '-r200')
    end
%     equalWeights = checkNURBSweightsCompatibility(nurbs);
end   
if prePlot.plot2Dgeometry
    figure('Color','white','name',['Cross section of Fluid 3D NURBS geometry. Mesh ' num2str(M)])
    for j = 1:numel(varCol)
        nurbs = subNURBS(varCol{j}.nurbs,'at',[1 0; 0 0; 0 0]);
        switch model
            case 'PH'
                prePlot.setViewAngle = [0,90];
            otherwise
                prePlot.setViewAngle = [0,0];
                switch varCol{j}.media
                    case 'fluid'
                        prePlot.color = [173, 216, 230]/255;
                    case 'solid'
                        prePlot.color = getColor(1);
                end
        end
        prePlot.displayName = ['Domain ' num2str(j)];
        plotNURBS(nurbs, prePlot);
    end
    hold off
    axis equal
    axis(prePlot.axis)
    if ~strcmp(prePlot.axis,'off')
        xlabel(prePlot.xlabel)
        ylabel(prePlot.ylabel)
        zlabel(prePlot.zlabel)
    end
    set(gca, 'Color', 'none');
    set(gca, 'Color', 'white');
    set(gcf,'color','w');
    drawnow
    if ~isempty(prePlot.export_fig_name2D)
        figureFullScreen(gcf)
        export_fig(['../../graphics/S1/SEM_p' num2str(degreeElev+4)], '-png', '-transparent', '-r200')
    end
            
    if ~isempty(prePlot.title)
        title(prePlot.title)
    end
    
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    if ~isempty(prePlot.export_fig_name2D)
        figureFullScreen(gcf)
        export_fig(prePlot.export_fig_name2D, '-png', '-transparent', '-r200')
    end
end



