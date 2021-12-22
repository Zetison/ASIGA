function task = plotMeshAndGeometry(task)
prePlot = task.prePlot;
if prePlot.plot3Dgeometry
    task.prePlot.fig3Dplot = figure('Color','white','name',['3D plot of geometry with mesh ' num2str(task.msh.M)]);
    if strcmp(task.misc.method, 'IENSG') && prePlot.plotArtificialBndry
        c_z = task.varCol{1}.c_z;
        c_xy = task.varCol{1}.c_xy;
        alignWithAxis = task.varCol{1}.alignWithAxis;
        x_0 = task.iem.x_0;
        ellipsoid = getEllipsoidData('C',[c_xy,c_xy,c_z]*task.iem.A_2,'alignWithAxis', alignWithAxis, 'x_0', x_0,'Xi',task.msh.Xi);
        if prePlot.alphaValue == 1
            prePlot.alphaValue = 0.8;
        end
        if prePlot.plotFullDomain
            plotNURBS(ellipsoid,'resolution',[100 100],'alphaValue',0.6,'color','blue');
        end
        noTopsets = numel(prePlot.plotSubsets);
        for i = 1:noTopsets
            if strcmp(prePlot.plotSubsets{i},'xy')
                theta = linspace(0,2*pi,10000);
                x = c_z*cos(theta);
                y = c_xy*sin(theta);
                plot(x,y,'black','DisplayName','Ellipse')
            end
        end
    end
    for j = 1:numel(task.varCol)
        if strcmp(task.misc.coreMethod, 'linear_FEM')
            prePlot.resolution = [0,0,0];
        end
        if isempty(prePlot.color)
            switch task.varCol{j}.media
                case 'fluid'
                    if task.varCol{j}.boundaryMethod
                        prePlot.color = getColor(1);
                    else
                        prePlot.color = getColor(10);
                    end
                case 'solid'
                    prePlot.color = getColor(1);
            end
        end
        nurbs = task.varCol{j}.nurbs;
        prePlot.displayName = ['Domain ' num2str(j)];
        if prePlot.plotFullDomain
            plotNURBS(nurbs, prePlot);
        end
        if isfield(task.varCol{j},'geometry') && ~isempty(prePlot.plotSubsets)
            topset = task.varCol{j}.geometry.topologysets.set;
            noTopsets = numel(prePlot.plotSubsets);
            colors = jet(noTopsets);
            colors(1:size(prePlot.color,1),:) = prePlot.color;
            for i = 1:noTopsets
                idx = findSet(topset,prePlot.plotSubsets{i});
                if isnan(idx)
                    warning(['Subset ' prePlot.plotSubsets{i} ' does not exist in structure'])
                    break
                end
                noPatches = numel(topset{idx}.item);
                nurbs = cell(1,noPatches);
                for ii = 1:numel(topset{idx}.item)
                    at = zeros(2,3,'logical');
                    patch = topset{idx}.item{ii}.Attributes.patch;
                    midx = topset{idx}.item{ii}.Text;
                    at(midx) = true;
                    outwardPointingNormals = strcmp(topset{idx}.Attributes.normal, 'outward');
                    inwardPointingNormals = strcmp(topset{idx}.Attributes.normal, 'inward');
    
                    nurbs(ii) = subNURBS(task.varCol{j}.nurbs(patch),'at',at.','outwardPointingNormals',outwardPointingNormals,'inwardPointingNormals',inwardPointingNormals);
                end
                prePlot.displayName = ['Domain ' num2str(j) ', ' topset{idx}.Attributes.name];
                if strcmp(topset{idx}.Attributes.name,'Gamma')
                    prePlot.color = getColor(1);
                else
                    prePlot.color = colors(i,:);
                end
                plotNURBS(nurbs, prePlot);
            end
        end
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
    figName = prePlot.export_fig_name3D;
    if isempty(figName)
        figName = task.saveName;
    end
    figName = [task.resultsFolder '/' figName '_3D'];
    if exist('../export_fig', 'dir')
        export_fig(figName, '-png', '-transparent', prePlot.pngResolution)
    end
    savefig([figName, '.fig'])
end   
if prePlot.plot2Dgeometry
    figure('Color','white','name',['Cross section of Fluid 3D NURBS geometry. Mesh ' num2str(task.msh.M)])
    for j = 1:numel(task.varCol)
        switch task.varCol{j}.media
            case 'fluid'
                prePlot.color = getColor(10);
            case 'solid'
                prePlot.color = getColor(1);
        end
        switch task.misc.model
            case {'MS'}
                nurbs2D = task.varCol{j}.nurbs(1:4:end);
                nurbs = subNURBS(nurbs2D,'at',[1 0; 0 0; 0 0]);
            otherwise
                nurbs2D = task.varCol{j}.nurbs;
                nurbs = subNURBS(nurbs2D,'at',[1 0; 0 0; 0 0]);
        end
        switch task.misc.model
            case {'PH', 'M3', 'MS', 'IMS', 'SMS','Mi2021ilc'}
                prePlot.view = [0,90];
            otherwise
                prePlot.view = [0,0];
        end
        if numel(task.varCol) > 1
            switch task.varCol{j}.media
                case 'fluid'
%                     prePlot.color = [1,1,1];
                    prePlot.color = getColor(10);
                case 'solid'
                    prePlot.color = getColor(1);
            end
        else
%             prePlot.color = [1,1,1];
            prePlot.color = getColor(10);
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
%     set(gca, 'Color', 'none');
    set(gca, 'Color', 'white');
    set(gcf,'color','w');
    drawnow
            
    if ~isempty(prePlot.title)
        title(prePlot.title)
    end
    
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    figName = prePlot.export_fig_name2D;
    if isempty(figName)
        figName = task.saveName;
    end
    figName = [task.resultsFolder '/' figName '_2D'];
    if exist('../export_fig', 'dir')
        ax.SortMethod='ChildOrder';
%         export_fig(figName, '-png', '-transparent', prePlot.pngResolution)
        export_fig(figName, '-pdf', '-transparent')
    end
    savefig([figName, '.fig'])
end



% plotControlPts(fluid);
% dofsToRemove = task.varCol{1}.dofsToRemove;
% controlPts = task.varCol{1}.controlPts;
% plot3(controlPts(dofsToRemove,1),controlPts(dofsToRemove,2),controlPts(dofsToRemove,3),'o','color','blue','MarkerFaceColor','blue')

